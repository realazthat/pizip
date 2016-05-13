#include <cstdint>
#include <vector>
#include <deque>
#include <string>
#include <ctime>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <bitset>
#include <thread>
#include <atomic>
#include <algorithm>
#include <cstring>
#include <cassert>

///threadpool lib (https://github.com/vit-vit/CTPL)
// #include "ctpl_stl.h"

///threadpool lib (https://github.com/maginatics/threadpool)
#include <maginatics/threadpool/threadpool.h>


std::size_t find_match_in_pi(const uint8_t* chunk_data, std::size_t chunk_size, const std::vector<uint8_t>& pi_data)
{
    for (std::size_t pi_index = 0; pi_index + chunk_size < pi_data.size(); ++pi_index)
    {

        if (memcmp(chunk_data, &pi_data[pi_index], chunk_size) == 0)
        {
            return pi_index;
        }
    }

    return pi_data.size();
}

template<typename T>
const char* print_time(std::ostream& out, std::clock_t& c_start, const T& t_start)
{
    std::clock_t c_end = std::clock();
    auto t_end = std::chrono::high_resolution_clock::now();

    out << std::fixed << std::setprecision(2) << "CPU: "
              << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << " ms, "
              << "Wall: "
              << std::chrono::duration<double, std::milli>(t_end-t_start).count()
              << " ms";

    return "";
}

void compress(    std::size_t threads
                , const uint8_t* unpacked_data, std::size_t unpacked_size
                , std::size_t chunk_size
                , std::ostream& out, const std::vector<uint8_t>& pi_data, std::ostream* log=nullptr)
{

    std::clock_t c_start = std::clock();
    auto t_start = std::chrono::high_resolution_clock::now();



    ///"0. Write the original input size as a 64-bit number"
    {
        uint64_t unpacked_size_ = unpacked_size;
        out.write(reinterpret_cast<const char*>(&unpacked_size_), sizeof(uint64_t));
    }

    if (log) *log << print_time(*log, c_start, t_start) << ", unpacked_size: " << unpacked_size << ", chunk_size: " << chunk_size << std::endl;

    ///"1) Divide the file into 16-byte chunks."
    std::size_t chunk_count = (unpacked_size / chunk_size) + ((unpacked_size % chunk_size) > 0 ? 1 : 0);


    ///every 8 chunks will have a byte (8 bits) to indicate if it has been compressed, or has not been compressed.
    ///thus the bitmap will have a size of approximately `chunk_count / 8`
    std::size_t bitmap_bytes = (chunk_count / 8) + ((chunk_count % 8) > 0 ? 1 : 0);

    if (log) *log << "unpacked_size: " << unpacked_size
                  << ", chunk_count: " << chunk_count
                  << ", bitmap_bytes: " << bitmap_bytes
                  << std::endl;

    ///"2) Create a bitmap with 1 bit per 16-byte chunk."
    typedef std::vector<uint8_t> bitmap_t;
    bitmap_t bitmap(bitmap_bytes,0);

    ///keep track of all the matches, and to where they match in the `pi_data`
    ///initialize the matches vector to all invalid match indices (`pi_data.size()` is an invalid match index)
    std::size_t invalid_match_index = pi_data.size();
    typedef std::deque<uint32_t> matches_t;
    matches_t matches(chunk_count, invalid_match_index);

    std::atomic<std::size_t> processed_count{0};
    std::atomic<std::size_t> matches_count{0};


    auto compute_chunk = [chunk_size, unpacked_data, unpacked_size, &pi_data, &processed_count, &matches_count, &matches](std::size_t chunk_index){
        const uint8_t* chunk_data = unpacked_data + chunk_index*chunk_size;
        std::size_t this_chunk_size = std::min(chunk_size, unpacked_size - chunk_size*chunk_index);

        assert(chunk_data+this_chunk_size <= unpacked_data + unpacked_size);
        std::size_t match_index = find_match_in_pi(chunk_data, this_chunk_size, pi_data);


        processed_count += 1;

        if (match_index < pi_data.size())
        {

            ///FOUND A MATCH!!

            matches_count += 1;


            
        }

        ///make sure pi_index fits into a 32 bit number
        assert (uint32_t(match_index) == match_index);

        matches.at(chunk_index) = match_index;
    };


    // ctpl::thread_pool p(threads);
    maginatics::ThreadPool p(threads);

    if (log) *log << print_time(*log, c_start, t_start) << ", starting threadpool with " << threads << " threads" << std::endl;

    auto print_stats = [&processed_count, &matches_count, chunk_count, &c_start, &t_start, log, &p](){

        if (!log)
            return;

        std::clock_t c_end = std::clock();
        auto t_end = std::chrono::high_resolution_clock::now();

        double c_time = 1000.0 * double(c_end-c_start) / CLOCKS_PER_SEC;
        double t_time = std::chrono::duration<double, std::milli>(t_end-t_start).count();


        double c_chunks_per_sec = 1000.0*double(processed_count) / c_time;
        double t_chunks_per_sec = 1000.0*double(processed_count) / t_time;
        double c_matches_per_sec = 1000.0*double(matches_count) / c_time;
        double t_matches_per_sec = 1000.0*double(matches_count) / t_time;
        *log << print_time(*log, c_start, t_start) << ", chunks processed: " << processed_count << "/" << chunk_count
                                           << ", matches_count: " << matches_count
                                           << ", CPU chunks/s: " << c_chunks_per_sec
                                           << ", chunks/s: " << t_chunks_per_sec
                                           << ", CPU matches/s: " << c_matches_per_sec
                                           << ", matches/s: " << t_matches_per_sec << std::endl
                                           << ", pool.queueLength(): " << p.queueLength() << std::endl;
    };



    for (std::size_t chunk_index = 0; chunk_index < chunk_count; ++chunk_index)
    {
        p.execute(std::bind(compute_chunk, chunk_index));
        print_stats();
    }

    while (p.queueLength() > 0)
    {
        print_stats();
        std::this_thread::sleep_for (std::chrono::seconds(3));
    }
    // p.join();

    for (std::size_t chunk_index = 0; chunk_index < chunk_count; ++chunk_index)
    {
        uint32_t match_index = matches.at(chunk_index);
        assert(match_index <= pi_data.size());
        if (!(match_index < pi_data.size()))
            continue;

        ///set the corresponding bit in the bitmap
        bitmap.at(chunk_index/8) |= (1 << (chunk_index % 8));
    }

    ///only keep the valid matches in the order they were found.
    auto last = std::remove(matches.begin(), matches.end(), invalid_match_index);

    matches.erase(last,matches.end());

    assert(matches.size() == matches_count);

    if (log) *log << print_time(*log, c_start, t_start) << ", matches.size(): " << matches.size() << std::endl;
    if (log) *log << print_time(*log, c_start, t_start) << ", now writing data out" << std::endl;

    ///now write the bitmap
    out.write(reinterpret_cast<const char*>(bitmap.data()), bitmap.size());

    ///now write out all the chunks
    for (std::size_t chunk_index = 0; chunk_index < chunk_count; ++chunk_index)
    {
        const uint8_t* chunk_data = unpacked_data + (chunk_index*chunk_size);
        
        ///the last chunk might be smaller
        std::size_t this_chunk_size = std::min(chunk_size, unpacked_size - chunk_size*chunk_index);

        bool has_match = bitmap.at(chunk_index/8) & (1 << (chunk_index % 8));

        if (has_match)
        {
            ///has match,
            ///"3) If the 16-byte sequence can be found in the first 2**32 digits of PI (encoding scheme TBD),
            /// set the bitmap index to 1 and write the 32-bit index."
            

            ///sanity: make sure there is a corresponding match left in `matches`
            assert(matches.size() > 0);

            ///retrieve the match index
            uint32_t match_index = matches.front(); matches.pop_front();

            assert(match_index < pi_data.size());

            ///has match,
            ///write the match index
            out.write(reinterpret_cast<const char*>(&match_index), sizeof(uint32_t));
        } else {
            ///no match,
            ///"4) If the 16-byte sequence cannot be found in the first 2**32 digits of PI (encoding scheme TBD),
            /// set the bitmap index to 0 and write the 16-byte chunk"
            
            ///write the chunk data
            ///sanity: make sure that we are not going to write from data outside the buffer.
            assert(chunk_data+this_chunk_size <= unpacked_data + unpacked_size);
            out.write(reinterpret_cast<const char*>(chunk_data), this_chunk_size);
        }
    }
    
    ///we should have consumed all of our matches while writing the chunks.
    assert(matches.size() == 0);

}

void uncompress( const uint8_t* packed_data, std::size_t packed_size, std::size_t chunk_size
                , std::ostream& out, const std::vector<uint8_t>& pi_data, std::ostream* log=nullptr)
{

    const uint8_t* cur = packed_data;

    if (!(cur + sizeof(uint64_t) < packed_data + packed_size))
    {
        std::cerr << "Unexpected end of data while reading UNPACKEDSIZE" << std::endl;
        exit(-1);
    }

    ///"0. Read the original input size as a 64-bit number"
    uint64_t unpacked_size = *reinterpret_cast<const uint64_t*>(cur);
    cur += sizeof(uint64_t);


    ///"1) Divide the file into 16-byte chunks."
    std::size_t chunk_count = (unpacked_size / chunk_size) + ((unpacked_size % chunk_count) > 0 ? 1 : 0);


    ///every 8 chunks will have a byte (8 bits) to indicate if it has been compressed, or has not been compressed.
    ///thus the bitmap will have a size of approximately `chunk_count / 8`
    std::size_t bitmap_bytes = chunk_count / 8 + ((chunk_count % 8) > 0 ? 1 : 0);

    if (log) *log << "cur: " << (cur-packed_data) << "/" << packed_size
                  << ", unpacked_size: " << unpacked_size
                  << ", chunk_count: " << chunk_count
                  << ", bitmap_bytes: " << bitmap_bytes
                  << std::endl;



    if (!(cur + bitmap_bytes < packed_data + packed_size))
    {
        std::cerr << "Unexpected end of data while reading BITMAP" << std::endl;
        exit(-1);
    }

    ///"2) Create a bitmap with 1 bit per 16-byte chunk."
    std::vector<uint8_t> bitmap(cur, cur+bitmap_bytes);
    cur += bitmap_bytes;

    auto count_bits = [](const std::vector<uint8_t>& bitmap){
        std::size_t onbits = 0;
        for (auto b : bitmap)
        {
            onbits += std::bitset<8>(b).count();
        }
        return onbits;
    };

    if (log) *log << "cur: " << (cur-packed_data) << "/" << packed_size
                  << ", bitmap matches: " << count_bits(bitmap)
                  << std::endl;

    ///now read in all the chunks
    for (std::size_t chunk_index = 0; chunk_index < chunk_count; ++chunk_index)
    {
        std::size_t this_chunk_size = std::min(chunk_size, unpacked_size - chunk_size*chunk_index);

        ///did this chunk match a part of the pi-sequence?
        bool had_match = bitmap.at(chunk_index/8) & (1 << (chunk_index % 8));

        assert(had_match == (bitmap.at(chunk_index/8) >> (chunk_index % 8)) & 1);

        if (log) *log << "cur: " << (cur-packed_data) << "/" << packed_size
                      << ", chunk_index: " << chunk_index
                      << ", this_chunk_size: " << this_chunk_size
                      << ", had_match: " << (had_match ? "true" : "false")
                      << std::endl;
        if (had_match)
        {

            if (!(cur + sizeof(uint32_t) < packed_data + packed_size))
            {
                std::cerr << "Unexpected end of data while reading MATCHINDEX" << std::endl;
                exit(-1);
            }

            uint32_t match_index = *reinterpret_cast<const uint32_t*>(cur);
            cur += sizeof(uint32_t);

            if (!(match_index + this_chunk_size <= pi_data.size()))
            {
                std::cerr << "Unexpected end of pi digit sequence; MATCHINDEX was: " << match_index
                          << ", with a chunk size of " << this_chunk_size
                           << ", but pi sequence only had " << pi_data.size() << " digits/bytes" << std::endl;
                exit(-1);
            }

            out.write(reinterpret_cast<const char*>(pi_data.data() + match_index), this_chunk_size);
        } else {

            if (!(cur + this_chunk_size < packed_data + packed_size))
            {
                std::cerr << "Unexpected end of data while reading chunk data" << std::endl;
                exit(-1);
            }
            out.write(reinterpret_cast<const char*>(cur), this_chunk_size);
            cur += this_chunk_size;
        }
    }

    ///sanity check: make sure the cur pointer reached the end of the incoming packed data.
    assert(cur == packed_data + packed_size);
}


void print_usage(std::ostream& out)
{
    out  << std::endl << "Usage: " << std::endl << std::endl
        << "  pizip [-d|--decompress] [-t <threads>]  [-c|--chunk <chunk-sizes>] </path/to/pihexdata> </path/to/infile> </path/to/outfile>" << std::endl
        << "    Compress or decompress a file (default is to compress)" << std::endl << std::endl
        << "      -d|--decompress            decompress the input file (default is to compress)" << std::endl << std::endl
        << "      -t <threads>               number of threads to use for compression (defaults to 1)" << std::endl << std::endl
        << "      -c|--chunk <chunk-sizes>   number of bytes per chunk (defaults to 16)" << std::endl << std::endl
        << "      /path/to/pihexdata         path to a file containing pi digits in hex" << std::endl << std::endl
        << "      /path/to/infile            path to input file" << std::endl << std::endl
        << "      /path/to/outfile           path to output file" << std::endl << std::endl
        << "  pizip [-h|--help]" << std::endl
        << "    show this help menu" << std::endl
        << std::endl
        << std::endl;
}

int main(int argc, const char** argv){

    if (argc == 2 && (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0))
    {
        print_usage(std::cerr);
        exit(0);
    }

    bool seen_decompress = false;
    bool decompress = false;

    bool seen_threads = false;
    std::size_t threads = 1;
    bool seen_chunk_size = false;
    std::size_t chunk_size = 16;

    const char* pihexdatapath = NULL;
    const char* infilepath = NULL;
    const char* outfilepath = NULL;

    ///parse the args
    for (std::size_t arg_index = 1; arg_index < argc; ++arg_index)
    {
        const char* arg = argv[arg_index];
        if (strcmp(arg, "-d") == 0 || strcmp(arg, "--decompress") == 0)
        {
            if (seen_decompress)
            {
                std::cerr << "Error: decompress can only be used once" << std::endl << std::endl;
                print_usage(std::cerr);
                exit(-1);
            }
            decompress = true;
            seen_decompress = true;
            continue;
        }
        if (strcmp(arg, "-t") == 0 || strcmp(arg, "--threads") == 0)
        {
            if (!(arg_index + 1 < argc))
            {
                std::cerr << "Error: expected number of threads; unexpected end of command" << std::endl << std::endl;
                print_usage(std::cerr);
                exit(-1);
            }

            if (seen_threads)
            {
                std::cerr << "Error: threads can only be specified once" << std::endl << std::endl;
                print_usage(std::cerr);
                exit(-1);
            }
            
            threads = strtoul(argv[arg_index + 1], NULL, 10);
            seen_threads = true;

            arg_index++;
            continue;
        }
        if (strcmp(arg, "-c") == 0 || strcmp(arg, "--chunk") == 0)
        {
            if (!(arg_index + 1 < argc))
            {
                std::cerr << "Error: expected chunk size; unexpected end of command" << std::endl << std::endl;
                print_usage(std::cerr);
                exit(-1);
            }

            if (seen_chunk_size)
            {
                std::cerr << "Error: chunk size can only be specified once" << std::endl << std::endl;
                print_usage(std::cerr);
                exit(-1);
            }
            

            chunk_size = strtoul(argv[arg_index + 1], NULL, 10);
            seen_chunk_size = true;

            arg_index++;
            continue;
        }

        if (!pihexdatapath)
        {
            pihexdatapath = arg;
            continue;
        }

        if (!infilepath)
        {
            infilepath = arg;
            continue;
        }
        if (!outfilepath)
        {
            outfilepath = arg;
            continue;
        }

        std::cerr << "Error: Too many arguments: \"" << arg << "\"" << std::endl << std::endl;
        print_usage(std::cerr);
        exit(-1);
    }


    std::ifstream pihexfile(pihexdatapath, std::ios::in | std::ios::binary);
    std::ifstream infile(infilepath, std::ios::in | std::ios::binary);
    std::ofstream outfile(outfilepath, std::ios::out | std::ios::binary);

    if (!pihexfile)
    {
        std::cerr << "Error opening pi hex digits file (path " << pihexdatapath << ")" << std::endl;
        exit(-1);
    }

    if (!infile)
    {
        std::cerr << "Error opening input file (path " << infilepath << ")" << std::endl;
        exit(-1);
    }

    if (!outfile)
    {
        std::cerr << "Error opening output file (path " << outfilepath << ")" << std::endl;
        exit(-1);
    }

    infile.seekg(0, infile.end);

    std::size_t file_size = infile.tellg();
    infile.seekg(0, infile.beg);

    pihexfile.seekg(0, pihexfile.end);

    std::size_t pi_hex_data_file_size = pihexfile.tellg();
    pihexfile.seekg(0, pihexfile.beg);



    std::vector<uint8_t> data0; data0.reserve(file_size);

    ///temp buffer for reading
    std::vector<uint8_t> buf(1024,0);




    std::vector<int8_t> pi_hex_data; pi_hex_data.reserve(pi_hex_data_file_size);

    
    while (pihexfile)
    {
        pihexfile.read(reinterpret_cast<char*>(buf.data()), buf.size());
        std::size_t bytesread = pihexfile.gcount();

        pi_hex_data.insert(pi_hex_data.end(), buf.begin(), buf.begin() + bytesread);
    }

    if (pi_hex_data.size() != pi_hex_data_file_size)
    {
        std::cerr << "pi hex digits data size (" << pi_hex_data.size() << ") does not match file size (" << pi_hex_data_file_size << ")... wat" << std::endl;
        exit(-1);
    }


    while (infile)
    {
        infile.read(reinterpret_cast<char*>(buf.data()), buf.size());
        std::size_t bytesread = infile.gcount();

        data0.insert(data0.end(), buf.begin(), buf.begin() + bytesread);
    }

    if (data0.size() != file_size)
    {
        std::cerr << "Data size (" << data0.size() << ") does not match file size (" << file_size << ")... wat" << std::endl;
        exit(-1);
    }


    std::vector<uint8_t> pi_data;

    for (std::size_t pi_hex_digit_index = 0; pi_hex_digit_index + 1 < pi_hex_data.size(); pi_hex_digit_index += 2)
    {

        char left_c = pi_hex_data[pi_hex_digit_index];
        char right_c = pi_hex_data[pi_hex_digit_index + 1];

        auto digit_to_int = [](char c){
            if (c >= '0' && c <= '9')
                return int(c - '0');
            if (c >= 'a' && c <= 'z')
                return int(c - 'a');
            if (c >= 'A' && c <= 'Z')
                return int(c - 'A');
            
            return -1;
        };

        int bin_digit = digit_to_int(left_c)*16 + digit_to_int(right_c);

        if (! (bin_digit >= 0 && bin_digit <= 255) )
        {
            std::cerr << "Invalid hex digit(s) in pi hex digits file: " << pi_hex_data[pi_hex_digit_index] << ", " << pi_hex_data[pi_hex_digit_index+1] << std::endl;
            exit(-1);
        }

        pi_data.push_back(uint8_t(bin_digit));
    }


    if (!decompress)
    {
        compress(threads, data0.data(), data0.size(), chunk_size, outfile, pi_data, &std::cerr);
        outfile.flush();
        outfile.close();
    } else {

        uncompress(data0.data(), data0.size(), chunk_size, outfile, pi_data, &std::cerr);
        outfile.flush();
        outfile.close();
    }

    // std::cerr << "pi_data.size(): " << pi_data.size() << std::endl;

    return 0;
}








