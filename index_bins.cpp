#include "algorithms.hpp"
#include <iostream>
#include "Utils/kmer.h"
#include <limits>
#include <string>
#include <queue>
#include <functional>
#include <cstdio>
#include <chrono>
#include "defaultColumn.hpp"
#include <omp.h>
#include "kDataframes/kDataFrameSTL.hpp"
#include "parallel_hashmap/phmap_dump.h"
#include <parallel_hashmap/btree.h>
#include <glob.h>
#include <ctime>

typedef std::chrono::high_resolution_clock Time;
using namespace std::chrono;

using LEGENDS_MAP = phmap::parallel_flat_hash_map<uint64_t,
    std::vector<uint32_t>,
    std::hash<uint64_t>,
    std::equal_to<uint64_t>,
    std::allocator<std::pair<const uint64_t, vector<uint32_t>>>,
    6>; // 6 submaps because colors will grow

using LEGENDS_MAP_OLD = phmap::parallel_flat_hash_map<uint64_t, std::vector<uint32_t>>;

using std::string;
using std::vector;
using std::cerr;
using std::cout;

using phmap::flat_hash_map;


inline std::vector<std::string> glob2(const std::string& pattern) {
    using namespace std;

    glob_t glob_result;
    memset(&glob_result, 0, sizeof(glob_result));

    int return_value = glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
    if (return_value != 0) {
        globfree(&glob_result);
        stringstream ss;
        ss << "glob() failed with return_value " << return_value << endl;
        throw std::runtime_error(ss.str());
    }

    vector<string> filenames;
    for (size_t i = 0; i < glob_result.gl_pathc; ++i)
        filenames.push_back(string(glob_result.gl_pathv[i]));
    globfree(&glob_result);
    return filenames;
}

void index_bins(string bins_dir, int kSize, kDataFrame* frame) {



    auto* colors = new StringColorColumn();
    frame->addColumn("color", (Column*)colors);
    flat_hash_map<string, string> namesMap;
    flat_hash_map<string, uint64_t> tagsMap;
    flat_hash_map<string, uint64_t> groupNameMap;
    auto* legend = new flat_hash_map<uint64_t, std::vector<uint32_t>>();
    flat_hash_map<uint64_t, uint64_t> colorsCount;
    uint64_t readID = 0, groupID = 1;
    priority_queue<uint64_t, vector<uint64_t>, std::greater<uint64_t>> freeColors;
    flat_hash_map<string, uint64_t> groupCounter;



    string seqName, groupName;
    int total_bins_number = 0;
    flat_hash_map<string, string> basename_to_path;
    for (const auto& dirEntry : glob2(bins_dir + "/*")) {
        string file_name = (string)dirEntry;
        size_t lastindex = file_name.find_last_of(".");
        string bin_prefix = file_name.substr(0, lastindex);
        std::string bin_basename = bin_prefix.substr(bin_prefix.find_last_of("/\\") + 1);


        std::string::size_type idx;
        idx = file_name.rfind('.');
        std::string extension = "";
        if (idx != std::string::npos) extension = file_name.substr(idx + 1);
        if (extension != "bin") {
            cerr << "skipping " << file_name << " does not have extension .bin" << endl;
            continue;
        }

        basename_to_path.insert(pair(bin_basename, file_name));

        total_bins_number++;

        seqName = bin_basename;
        groupName = bin_basename;

        namesMap.insert(make_pair(seqName, groupName));
        auto it = groupNameMap.find(groupName);
        groupCounter[groupName]++;
        if (it == groupNameMap.end()) {
            groupNameMap.insert(make_pair(groupName, groupID));
            tagsMap.insert(make_pair(to_string(groupID), groupID));
            vector<uint32_t> tmp;
            tmp.clear();
            tmp.push_back(groupID);
            legend->insert(make_pair(groupID, tmp));
            colorsCount.insert(make_pair(groupID, 0));
            groupID++;
        }
    }


    flat_hash_map<uint64_t, string> inv_groupNameMap;
    for (auto& _ : groupNameMap)
        inv_groupNameMap[_.second] = _.first;

    string kmer;

    readID = 0;
    int processed_bins_count = 0;
    auto begin_time = Time::now();
    uint_fast64_t current_kmers_numbers = 0;

    for (const auto& [bin_basename, bin_path] : basename_to_path) {
        cout << "Processing " << ++processed_bins_count << "/" << total_bins_number << " | " << bin_basename << " ... " << endl;

        flat_hash_map<uint64_t, uint64_t> convertMap;


        string readName = bin_basename;
        string groupName = bin_basename;

        auto it = namesMap.find(readName);
        if (it == namesMap.end()) {
            cerr << "WARNING: " << "read " << readName << "doesn't have a group. Please, check the names file." << endl;
            continue;
        }

        uint64_t readTag = groupNameMap.find(groupName)->second;


        convertMap.clear();
        convertMap.insert(make_pair(0, readTag));
        convertMap.insert(make_pair(readTag, readTag));

        begin_time = Time::now();
        phmap::flat_hash_set<uint64_t> bin_hashes;
        phmap::BinaryInputArchive ar_in(bin_path.c_str());
        bin_hashes.phmap_load(ar_in);
        for (const uint64_t& hashed_kmer : bin_hashes) {
            frame->insert(hashed_kmer);
            uint64_t kmerOrder = frame->getkmerOrder(hashed_kmer);
            if (colors->size() < frame->size() + 1)
                colors->resize(frame->size() + 1);
            uint64_t currentTag = colors->colors->index[kmerOrder];
            auto itc = convertMap.find(currentTag);
            if (itc == convertMap.end()) {
                vector<uint32_t> colors = legend->find(currentTag)->second;
                auto tmpiT = find(colors.begin(), colors.end(), readTag);
                if (tmpiT == colors.end()) {
                    colors.push_back(readTag);
                    sort(colors.begin(), colors.end());
                }

                string colorsString = to_string(colors[0]);
                for (unsigned int k = 1; k < colors.size(); k++) {
                    colorsString += ";" + to_string(colors[k]);
                }

                auto itTag = tagsMap.find(colorsString);
                if (itTag == tagsMap.end()) {
                    uint64_t newColor;
                    if (freeColors.empty()) {
                        newColor = groupID++;
                    }
                    else {
                        newColor = freeColors.top();
                        freeColors.pop();
                    }

                    tagsMap.insert(make_pair(colorsString, newColor));
                    legend->insert(make_pair(newColor, colors));
                    itTag = tagsMap.find(colorsString);
                    colorsCount[newColor] = 0;
                    // if(groupID>=maxTagValue){
                    //   cerr<<"Tag size is not enough. ids reached "<<groupID<<endl;
                    //   return -1;
                    // }
                }
                uint64_t newColor = itTag->second;

                convertMap.insert(make_pair(currentTag, newColor));
                itc = convertMap.find(currentTag);
            }

            if (itc->second != currentTag) {

                colorsCount[currentTag]--;
                if (colorsCount[currentTag] == 0 && currentTag != 0) {
                    auto _invGroupNameIT = inv_groupNameMap.find(currentTag);
                    if (_invGroupNameIT == inv_groupNameMap.end()) {
                        freeColors.push(currentTag);

                        vector<uint32_t> colors = legend->find(currentTag)->second;
                        string colorsString = to_string(colors[0]);
                        for (unsigned int k = 1; k < colors.size(); k++) {
                            colorsString += ";" + to_string(colors[k]);
                        }
                        tagsMap.erase(colorsString);

                        legend->erase(currentTag);
                        if (convertMap.find(currentTag) != convertMap.end())
                            convertMap.erase(currentTag);
                    }

                }
                colorsCount[itc->second]++;
            }

            colors->colors->index[kmerOrder] = itc->second;

        }
        readID += 1;
        groupCounter[groupName]--;
        if (colorsCount[readTag] == 0) {
            if (groupCounter[groupName] == 0) {
                freeColors.push(readTag);
                legend->erase(readTag);
            }
        }

        auto loop_time_secs = std::chrono::duration<double, std::milli>(Time::now() - begin_time).count() / 1000;
        cout << "   loaded_kmers      " << bin_hashes.size() << endl;
        cout << "   uniq_added_kmers: " << frame->size() - current_kmers_numbers << endl;
        cout << "   total_kmers       " << frame->size() << " | load_factor: " << frame->load_factor() << endl;
        cout << "   total_colors      " << legend->size() << endl;
        cout << "   loop_time:        " << loop_time_secs << " secs" << endl;
        cout << "--------" << endl;
        current_kmers_numbers = frame->size();

    }

    colors->colors->values = new mixVectors(*legend, groupCounter.size());
    delete legend;
    for (auto& iit : namesMap) {
        uint32_t sampleID = groupNameMap[iit.second];
        colors->namesMap[sampleID] = iit.second;
    }
    frame->colorColumn = (colors->colors);

}

inline uint64_t to_uint64_t(std::string const& value) {
  uint64_t result = 0;
  char const* p = value.c_str();
  char const* q = p + value.size();
  while (p < q) {
    result *= 10;
    result += *(p++) - '0';
  }
  return result;
}

int main(int argc, char** argv) {
    if(argc < 5){
        throw "args: <bins_dir> <kSize> <output_prefix> <initial_reserve_size>\n";
    }
    string bins_dir = argv[1];
    int kSize = stoi(argv[2]);
    string output_prefix = argv[3];
    uint64_t reserve_size = to_uint64_t(argv[4]);


    auto* kf = kDataFrameFactory::createPHMAP(kSize, reserve_size);
    index_bins(bins_dir, kSize, kf);
    kf->save(output_prefix);
}