// Streaming, bounded-memory screen for sequence-visible PacBio CCS artefacts.
// C++17; no third-party libraries. Decisions are automated FASTA heuristics.
// Developed by GPT-6 Sol under the guidance and supervision of Yanpeng Chen.
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cctype>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <limits>
#include <exception>
#include <mutex>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>
#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

struct Record { std::string id, header, seq; };
struct Adapter { std::string name, seq, strand; };
struct Options {
    std::string input, output, summary, adapters, retained_fasta;
    unsigned log_every = 5;
    unsigned threads = std::min(8u, std::max(1u, std::thread::hardware_concurrency()));
    size_t batch_bases = 2 * 1024 * 1024;
    int min_repeat_span = 1800;
    int min_seed_pairs = 8;
    int seed_mod = 16;
    bool all = false, force = false;
};
struct Counters {
    uint64_t reads = 0, bases = 0, exclude_adapter = 0, exclude_palindrome = 0;
    uint64_t quarantine = 0, keep_direct_repeat = 0, keep = 0;
    uint64_t adapter = 0, direct = 0, inverted = 0;
    Counters& operator+=(const Counters& b) {
        reads += b.reads; bases += b.bases; exclude_adapter += b.exclude_adapter;
        exclude_palindrome += b.exclude_palindrome; quarantine += b.quarantine;
        keep_direct_repeat += b.keep_direct_repeat; keep += b.keep;
        adapter += b.adapter; direct += b.direct;
        inverted += b.inverted; return *this;
    }
};
struct Batch { uint64_t index; std::vector<Record> reads; };
struct BatchResult { uint64_t index; std::string rows, quarantine_ids, direct_repeat_ids, retained_fasta; Counters counts; };

static void usage() {
#ifdef _WIN32
    const bool tty = _isatty(_fileno(stdout));
#else
    const bool tty = isatty(fileno(stdout));
#endif
    const char* term = std::getenv("TERM");
    const bool color = tty && !std::getenv("NO_COLOR") && (!term || std::string(term)!="dumb");
    const char* heading = color ? "\033[1;36m" : "";
    const char* reset = color ? "\033[0m" : "";
    std::cout << heading << "HiFiCheck | streaming FASTA artefact screen" << reset << "\n"
      "Usage: hificheck --input reads.fa --output results.tsv [options]\n\n"
      << heading << "Input and output" << reset << "\n"
      "  -i, --input FILE       Plain FASTA file; '-' reads stdin. No direct gzip support.\n"
      "  -o, --output FILE      Exact TSV file path, NOT a directory or filename prefix.\n"
      "                        Only flagged reads by default; --all includes every read.\n"
      "  --summary FILE        Exact JSON file path (default: OUTPUT.summary.json)\n"
      "  --retained-fasta FILE  Exact FASTA file path for keep + keep_direct_repeat reads\n"
      "  --adapters FILE        Library-specific adapter FASTA; optional\n"
      "  --all                  Include unflagged reads in TSV\n"
      "  --force                Replace existing output files for a new run\n\n"
      << heading << "Performance and logging" << reset << "\n"
      "  -t, --threads N        Worker threads (default: min(8, hardware concurrency))\n"
      "  --batch-bases N        Approx. bases per work batch (default: 2097152)\n"
      "  --log-every N          Progress interval in seconds (default: 5; 0 disables)\n"
      "                        Logs go to stderr; counts, elapsed time and throughput.\n\n"
      << heading << "Screening parameters" << reset << "\n"
      "  --min-repeat-span N    Seeded repeat span in bp (default: 1800)\n"
      "  --min-seed-pairs N     Supporting seed pairs (default: 8)\n"
      "  --seed-mod N           Subsample 17-mers by factor N, power of 2 (default: 16)\n"
      "  -h, --help             Show this help (terminal colors; NO_COLOR disables)\n\n"
      << heading << "Output example: --output sample.hificheck.tsv" << reset << "\n"
      "  sample.hificheck.tsv                  Read-level decisions\n"
      "  sample.hificheck.tsv.summary.json     Counts, parameters and throughput\n"
      "  sample.hificheck.tsv.quarantine.ids   Excluded + quarantined read IDs\n"
      "  sample.hificheck.tsv.direct_repeat.ids  Retained direct-repeat candidate IDs\n"
      "Suffixes are appended to the ENTIRE output path, including .tsv.\n"
      "--summary overrides the JSON path; --retained-fasta creates an optional file.\n\n"
      << heading << "Gzip input (Bash)" << reset << "\n"
      "  set -o pipefail\n"
      "  gzip -dc reads.fa.gz | hificheck -i - -o results.tsv -t 24\n\n"
      "Coordinates are zero-based, half-open. Adapter matching permits substitutions\n"
      "up to 10% of adapter length; indels and BAM/quality evidence are not assessed.\n";
}
static unsigned long long parse_ull(const std::string& value, const char* label) {
    if (value.empty() || value.find_first_not_of("0123456789") != std::string::npos)
        throw std::runtime_error(std::string("Invalid ") + label + ": " + value);
    size_t end = 0;
    unsigned long long n = std::stoull(value, &end);
    if (end != value.size()) throw std::runtime_error(std::string("Invalid ") + label + ": " + value);
    if (n > static_cast<unsigned long long>(std::numeric_limits<int>::max()))
        throw std::runtime_error(std::string("Numeric parameter out of range: ") + label);
    return n;
}
static bool same_path(const std::filesystem::path& a, const std::filesystem::path& b) {
    return std::filesystem::weakly_canonical(a)==std::filesystem::weakly_canonical(b) ||
        (std::filesystem::exists(a) && std::filesystem::exists(b) && std::filesystem::equivalent(a,b));
}
static Options parse_args(int argc, char** argv) {
    Options o;
    for (int i=1; i<argc; ++i) {
        std::string k = argv[i];
        if (k == "-h" || k == "--help") { usage(); std::exit(0); }
        if (k == "--all") { o.all=true; continue; }
        if (k == "--force") { o.force=true; continue; }
        if (i+1 >= argc) throw std::runtime_error("Missing value for " + k);
        std::string v = argv[++i];
        if (k == "-i" || k == "--input") o.input=v;
        else if (k == "-o" || k == "--output") o.output=v;
        else if (k == "--summary") o.summary=v;
        else if (k == "--adapters") o.adapters=v;
        else if (k == "--retained-fasta") o.retained_fasta=v;
        else if (k == "--log-every") o.log_every=static_cast<unsigned>(parse_ull(v,"log-every"));
        else if (k == "-t" || k == "--threads") o.threads=static_cast<unsigned>(parse_ull(v,"threads"));
        else if (k == "--min-repeat-span") o.min_repeat_span=static_cast<int>(parse_ull(v,"min-repeat-span"));
        else if (k == "--min-seed-pairs") o.min_seed_pairs=static_cast<int>(parse_ull(v,"min-seed-pairs"));
        else if (k == "--seed-mod") o.seed_mod=static_cast<int>(parse_ull(v,"seed-mod"));
        else if (k == "--batch-bases") o.batch_bases=static_cast<size_t>(parse_ull(v,"batch-bases"));
        else throw std::runtime_error("Unknown option: " + k);
    }
    if (o.input.empty() || o.output.empty()) throw std::runtime_error("--input and --output are required");
    if (o.input=="-" && o.adapters=="-") throw std::runtime_error("Input and adapters cannot both use stdin");
    if (o.summary.empty()) o.summary=o.output+".summary.json";
    std::vector<std::string> output_paths={o.output,o.summary,o.output+".quarantine.ids",
        o.output+".direct_repeat.ids"};
    if(!o.retained_fasta.empty()) output_paths.push_back(o.retained_fasta);
    for(size_t j=0;j<output_paths.size();++j) {
        const auto path=std::filesystem::absolute(output_paths[j]).lexically_normal();
        if(o.input!="-" && same_path(path,o.input))
            throw std::runtime_error("Output path overlaps input FASTA");
        if(!o.adapters.empty() && o.adapters!="-" && same_path(path,o.adapters))
            throw std::runtime_error("Output path overlaps adapter FASTA");
        for(size_t k=0;k<j;++k) if(same_path(path,output_paths[k]))
            throw std::runtime_error("Output paths must differ");
        if(!o.force && std::filesystem::exists(path))
            throw std::runtime_error("Output already exists (use --force to replace): "+output_paths[j]);
    }
    if (o.threads < 1 || o.threads > 256 || o.batch_bases < 1000 ||
        o.min_repeat_span < 100 || o.min_seed_pairs < 2 ||
        o.seed_mod < 1 || o.seed_mod > 256 || (o.seed_mod & (o.seed_mod-1)))
        throw std::runtime_error("Invalid numeric parameter (see --help)");
    return o;
}
static int base_code(char c) {
    switch(c) { case 'A': return 0; case 'C': return 1; case 'G': return 2;
                case 'T': return 3; default: return -1; }
}
static std::string reverse_complement(const std::string& s) {
    std::string r(s.rbegin(), s.rend());
    for (char& c:r) {
        switch(c) { case 'A': c='T'; break; case 'C': c='G'; break;
                    case 'G': c='C'; break; case 'T': c='A'; break; default: c='N'; }
    }
    return r;
}
static std::string clean_id(const std::string& h) {
    size_t end=h.find_first_of(" \t\r");
    if (h.empty() || end==0) throw std::runtime_error("Empty FASTA identifier");
    std::string id=h.substr(0,end);
    for(char c:id) if(c=='\t' || c=='\n' || c=='\r') throw std::runtime_error("Invalid FASTA identifier");
    return id;
}
static void append_sequence(std::string& seq, const std::string& line) {
    for(char c:line) {
        if(c=='\r') continue;
        if(c==' ' || c=='\t') throw std::runtime_error("Whitespace within FASTA sequence");
        c=static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        if(c!='A' && c!='C' && c!='G' && c!='T' && c!='N') {
            if(std::string("RYWSKMBDHV").find(c)==std::string::npos)
                throw std::runtime_error("FASTA must contain IUPAC DNA bases only");
            // Preserve IUPAC symbols; base_code rejects them as repeat seeds.
        }
        seq.push_back(c);
    }
}
class FastaReader {
    std::ifstream file_;
    std::istream* input_;
    std::string header_;
    bool have_header_=false;
public:
    explicit FastaReader(const std::string& path) : input_(&std::cin) {
        if (path!="-") {
            file_.open(path, std::ios::binary);
            if (!file_) throw std::runtime_error("Cannot open FASTA: " + path);
            input_=&file_;
        }
    }
    bool next(Record& r) {
        std::string line;
        if (!have_header_) {
            while(std::getline(*input_,line)) {
                if(line.empty()) continue;
                if(line[0]!='>') throw std::runtime_error("Expected FASTA header (>), got sequence/data");
                header_=line.substr(1); have_header_=true; break;
            }
            if(!have_header_) {
                if(input_->bad()) throw std::runtime_error("FASTA read error");
                return false;
            }
        }
        r.id=clean_id(header_); r.header=header_; r.seq.clear();
        while(std::getline(*input_,line)) {
            if(!line.empty() && line[0]=='>') {
                header_=line.substr(1);
                if(r.seq.empty()) throw std::runtime_error("Empty sequence: " + r.id);
                return true;
            }
            append_sequence(r.seq,line);
        }
        if(input_->bad()) throw std::runtime_error("FASTA read error");
        have_header_=false;
        if(r.seq.empty()) throw std::runtime_error("Empty sequence: " + r.id);
        return true;
    }
};
static std::vector<Adapter> load_adapters(const std::string& path) {
    std::vector<Adapter> a;
    if(path.empty()) return a;
    FastaReader reader(path); Record r;
    while(reader.next(r)) {
        if(r.seq.size()<24 || r.seq.find_first_not_of("ACGT")!=std::string::npos)
            throw std::runtime_error("Adapter must have >=24 A/C/G/T bases: " + r.id);
        a.push_back({r.id,r.seq,"+"});
        a.push_back({r.id,reverse_complement(r.seq),"-"});
    }
    if(a.empty()) throw std::runtime_error("Adapter FASTA is empty");
    return a;
}
struct AdapterHit { bool found=false, internal=false; std::string name,strand; size_t start=0,end=0; int mismatches=0; };
static AdapterHit scan_adapters(const std::string& seq, const std::vector<Adapter>& adapters) {
    AdapterHit best;
    for(const auto& a:adapters) {
        const size_t len=a.seq.size();
        if(seq.size()<len) continue;
        const int maxmm=static_cast<int>(len/10);
        const int blocks=maxmm+1; // Pigeonhole guarantee for substitution-only matches.
        for(int b=0;b<blocks;++b) {
            const size_t from=len*static_cast<size_t>(b)/blocks;
            const size_t to=len*static_cast<size_t>(b+1)/blocks;
            const std::string seed=a.seq.substr(from,to-from);
            size_t pos=seq.find(seed);
            while(pos!=std::string::npos) {
                if(pos>=from) {
                    const size_t start=pos-from;
                    if(start+len<=seq.size()) {
                        int mm=0;
                        for(size_t j=0;j<len && mm<=maxmm;++j) mm+=(seq[start+j]!=a.seq[j]);
                        if(mm<=maxmm) {
                            const bool internal=(start>=100 && start+len+100<=seq.size());
                            if(!best.found || (internal && !best.internal) ||
                               (internal==best.internal && mm<best.mismatches))
                                best={true,internal,a.name,a.strand,start,start+len,mm};
                            if(internal && mm==0) return best;
                        }
                    }
                }
                pos=seq.find(seed,pos+1);
            }
        }
    }
    return best;
}
static uint64_t mix64(uint64_t x) {
    x+=0x9e3779b97f4a7c15ULL;
    x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;
    x=(x^(x>>27))*0x94d049bb133111ebULL;
    return x^(x>>31);
}
struct SeedPos { uint32_t pos; bool reverse; };
struct Bin { int pairs=0; uint32_t lo=UINT32_MAX,hi=0, other_lo=UINT32_MAX,other_hi=0; };
struct RepeatHit { bool found=false; char orientation='.'; uint32_t first_start=0,first_end=0,second_start=0,second_end=0; int pairs=0; };
static RepeatHit scan_repeat(const std::string& seq, const Options& o,
     std::unordered_map<uint64_t,std::vector<SeedPos>>& index,
     std::unordered_map<uint64_t,Bin>& bins) {
    constexpr int K=17;
    if(seq.size()<static_cast<size_t>(o.min_repeat_span)*2) return {};
    index.clear(); bins.clear();
    const uint64_t mask=(1ULL<<(2*K))-1;
    uint64_t forward=0,reverse=0;
    int valid=0;
    if(seq.size()>UINT32_MAX/2) throw std::runtime_error("Read exceeds supported coordinate range");
    for(uint32_t i=0;i<seq.size();++i) {
        int b=base_code(seq[i]);
        if(b<0) { forward=reverse=0; valid=0; continue; }
        forward=((forward<<2)|static_cast<uint64_t>(b))&mask;
        reverse=(reverse>>2)|(static_cast<uint64_t>(3-b)<<(2*(K-1)));
        if(valid<K) ++valid;
        if(valid<K) continue;
        const uint64_t canonical=std::min(forward,reverse);
        if((mix64(canonical)&static_cast<uint64_t>(o.seed_mod-1))!=0) continue;
        const bool rev=(reverse<forward);
        const uint32_t pos=i-K+1;
        auto& prior=index[canonical];
        if(prior.size()<=16) {
            for(const auto& p:prior) {
                const bool inverted=(p.reverse!=rev);
                // Opposite-orientation arms may meet at a short adapter-sized
                // junction; imposing the direct-repeat spacing loses that arm.
                if(pos-p.pos<static_cast<uint32_t>(inverted ? K : o.min_repeat_span)) continue;
                const uint32_t distance=inverted ? pos+p.pos : pos-p.pos;
                const uint64_t key=(static_cast<uint64_t>(distance/100)<<1)|static_cast<uint64_t>(inverted);
                Bin& bin=bins[key];
                ++bin.pairs;
                bin.lo=std::min(bin.lo,p.pos); bin.hi=std::max(bin.hi,p.pos);
                bin.other_lo=std::min(bin.other_lo,pos); bin.other_hi=std::max(bin.other_hi,pos);
            }
        }
        if(prior.size()<17) prior.push_back({pos,rev});
    }
    RepeatHit best;
    for(const auto& kv:bins) {
        const Bin& b=kv.second;
        const int span=static_cast<int>(b.hi-b.lo+K);
        if(b.pairs<o.min_seed_pairs || span<o.min_repeat_span) continue;
        const char orientation=(kv.first&1)?'-':'+';
        if(!best.found || (orientation=='-' && best.orientation=='+') ||
           (orientation==best.orientation && (b.pairs>best.pairs ||
            (b.pairs==best.pairs && b.lo<best.first_start)))) {
            best={true,(kv.first&1)?'-':'+',b.lo,b.hi+K,b.other_lo,b.other_hi+K,b.pairs};
        }
    }
    return best;
}
static void add_row(std::string& rows,const Record& r,const AdapterHit& a,const RepeatHit& rep,const char* action) {
    rows+=r.id; rows+='\t'; rows+=std::to_string(r.seq.size()); rows+='\t'; rows+=action; rows+='\t';
    rows+=(a.found ? (a.internal ? "internal_adapter" : "terminal_adapter") : "");
    if(rep.found) { if(a.found) rows+='|'; rows+=(rep.orientation=='+' ? "direct_repeat" : "inverted_repeat"); }
    rows+='\t'; rows+=a.found?a.name:""; rows+='\t'; rows+=a.found?a.strand:""; rows+='\t';
    if(a.found) { rows+=std::to_string(a.start); rows+='\t'; rows+=std::to_string(a.end); rows+='\t'; rows+=std::to_string(a.mismatches); }
    else rows+="\t\t";
    rows+='\t'; if(rep.found) rows+=rep.orientation; rows+='\t';
    if(rep.found) {
        rows+=std::to_string(rep.first_start); rows+='\t'; rows+=std::to_string(rep.first_end); rows+='\t';
        rows+=std::to_string(rep.second_start); rows+='\t'; rows+=std::to_string(rep.second_end); rows+='\t';
        rows+=std::to_string(rep.pairs);
    } else rows+="\t\t\t\t";
    rows+='\n';
}
static BatchResult process_batch(Batch batch,const Options& o,const std::vector<Adapter>& adapters) {
    BatchResult out; out.index=batch.index;
    std::unordered_map<uint64_t,std::vector<SeedPos>> index;
    std::unordered_map<uint64_t,Bin> bins;
    for(const auto& r:batch.reads) {
        AdapterHit a=scan_adapters(r.seq,adapters);
        RepeatHit rep=scan_repeat(r.seq,o,index,bins);
        ++out.counts.reads; out.counts.bases+=r.seq.size();
        if(a.found) ++out.counts.adapter;
        if(rep.found) {
            if(rep.orientation=='+') ++out.counts.direct;
            else ++out.counts.inverted;
        }
        // PacBio's adapter-palindrome signature has a near-terminal reverse-
        // complement arm. Other inverted matches remain uncertain from FASTA.
        const uint64_t first_span=rep.found ? rep.first_end-rep.first_start : 0;
        const uint64_t junction_gap=rep.found && rep.second_start>=rep.first_end ?
            rep.second_start-rep.first_end : UINT64_MAX;
        const bool terminal_palindrome=rep.found && rep.orientation=='-' &&
            (rep.first_start<=200 || rep.second_end+200>=r.seq.size()) &&
            2*first_span>=0.7*r.seq.size() && junction_gap<=500;
        const char* action;
        if(a.internal) { action="exclude_adapter"; ++out.counts.exclude_adapter; }
        else if(terminal_palindrome) { action="exclude_palindrome"; ++out.counts.exclude_palindrome; }
        else if(a.found || (rep.found && rep.orientation=='-')) {
            action="quarantine"; ++out.counts.quarantine;
        } else if(rep.found) {
            action="keep_direct_repeat"; ++out.counts.keep_direct_repeat;
        } else { action="keep"; ++out.counts.keep; }
        const std::string_view action_view(action);
        if(action_view=="exclude_adapter" || action_view=="exclude_palindrome" || action_view=="quarantine") {
            out.quarantine_ids+=r.id; out.quarantine_ids+='\n';
        }
        if(action_view=="keep_direct_repeat") {
            out.direct_repeat_ids+=r.id; out.direct_repeat_ids+='\n';
        }
        if(!o.retained_fasta.empty() && (action_view=="keep" || action_view=="keep_direct_repeat")) {
            out.retained_fasta+='>'; out.retained_fasta+=r.header;
            out.retained_fasta+='\n'; out.retained_fasta+=r.seq; out.retained_fasta+='\n';
        }
        if(o.all || action_view!="keep") add_row(out.rows,r,a,rep,action);
    }
    return out;
}
static void json_string(std::ostream& o,const std::string& s) {
    o<<'"'; for(char c:s) { if(c=='\\'||c=='"') o<<'\\'; if(static_cast<unsigned char>(c)<32) o<<'?'; else o<<c; } o<<'"';
}
// Independent timer keeps reporting even while stdin or a worker is blocked.
// Counters distinguish parsed input from results actually written in order.
class Progress {
    std::mutex mutex_;
    std::condition_variable cv_;
    bool stopped_=false;
    std::thread timer_;
public:
    std::atomic<uint64_t> loaded{0}, loaded_bases{0}, processed{0}, processed_bases{0};
    std::atomic<bool> eof{false};
    explicit Progress(unsigned interval) {
        if(!interval) return;
        timer_=std::thread([this,interval] {
            const auto start=std::chrono::steady_clock::now();
            std::unique_lock<std::mutex> lock(mutex_);
            while(!cv_.wait_for(lock,std::chrono::seconds(interval),[this]{return stopped_;})) {
                const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
                std::cerr << "[progress] elapsed=" << std::fixed << std::setprecision(1) << seconds
                          << "s loaded=" << loaded.load() << " reads/" << loaded_bases.load()
                          << " bases processed=" << processed.load() << " reads/" << processed_bases.load()
                          << " bases rate=" << processed.load()/seconds << " reads/s input="
                          << (eof.load()?"EOF":"reading") << "\n";
            }
        });
    }
    void stop() {
        { std::lock_guard<std::mutex> lock(mutex_); stopped_=true; }
        cv_.notify_all();
        if(timer_.joinable()) timer_.join();
    }
    ~Progress() { stop(); }
};
int main(int argc,char** argv) {
    try {
        const Options o=parse_args(argc,argv);
        const auto started=std::chrono::steady_clock::now();
        // A forced new run removes an old success marker before touching output.
        if(o.force) std::filesystem::remove(o.summary);
        std::cerr << "[start] input=" << (o.input=="-"?"stdin (plain FASTA)":o.input)
                  << " threads=" << o.threads << " batch_bases=" << o.batch_bases
                  << " log_every=" << o.log_every << "s\n"
                  << "[output] TSV file: " << o.output << "\n"
                  << "[output] JSON summary: " << o.summary << "\n"
                  << "[output] Excluded/quarantined IDs: " << o.output << ".quarantine.ids\n"
                  << "[output] Retained direct-repeat IDs: " << o.output << ".direct_repeat.ids\n";
        if(!o.retained_fasta.empty()) std::cerr << "[output] Retained FASTA: " << o.retained_fasta << "\n";
        std::cerr << "[screen] adapters=" << (o.adapters.empty()?"unassessed (no reference supplied)":o.adapters)
                  << "; TSV rows=" << (o.all?"all reads":"flagged reads only") << "\n";
        const auto adapters=load_adapters(o.adapters);
        FastaReader reader(o.input);
        std::ofstream output(o.output,std::ios::binary);
        if(!output) throw std::runtime_error("Cannot open output: "+o.output);
        std::ofstream quarantine_ids(o.output+".quarantine.ids",std::ios::binary);
        std::ofstream direct_repeat_ids(o.output+".direct_repeat.ids",std::ios::binary);
        if(!quarantine_ids || !direct_repeat_ids) throw std::runtime_error("Cannot open ID outputs");
        std::ofstream retained_fasta;
        if(!o.retained_fasta.empty()) {
            retained_fasta.open(o.retained_fasta,std::ios::binary);
            if(!retained_fasta) throw std::runtime_error("Cannot open retained FASTA: "+o.retained_fasta);
        }
        output<<"read_id\tread_length\tqc_action\tflags\tadapter_name\tadapter_strand\tadapter_start_0\tadapter_end_0\tadapter_mismatches\trepeat_orientation\trepeat_first_start_0\trepeat_first_end_0\trepeat_second_start_0\trepeat_second_end_0\trepeat_seed_pairs\n";
        // A fixed worker pool and capped in-flight work keep memory bounded.
        std::mutex mutex;
        std::condition_variable task_cv,done_cv;
        std::deque<Batch> tasks;
        std::map<uint64_t,BatchResult> completed;
        bool done=false;
        std::exception_ptr worker_error;
        Progress progress(o.log_every);
        std::vector<std::thread> workers;
        try {
        for(unsigned t=0;t<o.threads;++t) workers.emplace_back([&]() {
            try {
            while(true) {
                Batch b;
                {
                    std::unique_lock<std::mutex> lock(mutex);
                    task_cv.wait(lock,[&]{return done||!tasks.empty();});
                    if(worker_error || tasks.empty()) return;
                    b=std::move(tasks.front()); tasks.pop_front();
                }
                BatchResult result=process_batch(std::move(b),o,adapters);
                {
                    std::lock_guard<std::mutex> lock(mutex);
                    completed.emplace(result.index,std::move(result));
                }
                done_cv.notify_one();
            }
            } catch(...) {
                { std::lock_guard<std::mutex> lock(mutex);
                  if(!worker_error) worker_error=std::current_exception();
                  done=true; tasks.clear(); }
                task_cv.notify_all(); done_cv.notify_all();
            }
        });
        } catch(...) {
            { std::lock_guard<std::mutex> lock(mutex); done=true; }
            task_cv.notify_all();
            for(auto& worker:workers) if(worker.joinable()) worker.join();
            throw;
        }
        Counters total;
        uint64_t issued=0, written=0;
        auto collect=[&]() {
            BatchResult result;
            {
                std::unique_lock<std::mutex> lock(mutex);
                done_cv.wait(lock,[&]{return worker_error || completed.count(written)!=0;});
                if(worker_error) std::rethrow_exception(worker_error);
                auto it=completed.find(written);
                result=std::move(it->second); completed.erase(it);
            }
            output<<result.rows;
            quarantine_ids<<result.quarantine_ids;
            direct_repeat_ids<<result.direct_repeat_ids;
            if(retained_fasta.is_open()) retained_fasta<<result.retained_fasta;
            if(!output || !quarantine_ids || !direct_repeat_ids ||
               (retained_fasta.is_open() && !retained_fasta)) throw std::runtime_error("Failed writing output");
            total+=result.counts; ++written;
            progress.processed=total.reads;
            progress.processed_bases=total.bases;
        };
        try {
            Batch batch; size_t batch_bases=0; Record r;
            const uint64_t max_inflight=std::max<uint64_t>(2,o.threads*3ULL);
            while(reader.next(r)) {
                ++progress.loaded; progress.loaded_bases+=r.seq.size();
                batch_bases+=r.seq.size();
                batch.reads.push_back(std::move(r)); r=Record{};
                if(batch_bases>=o.batch_bases) {
                    batch.index=issued++;
                    { std::lock_guard<std::mutex> lock(mutex); tasks.push_back(std::move(batch)); }
                    task_cv.notify_one(); batch=Batch{}; batch_bases=0;
                    if(issued-written>=max_inflight) collect();
                }
            }
            progress.eof=true;
            if(!batch.reads.empty()) {
                batch.index=issued++;
                { std::lock_guard<std::mutex> lock(mutex); tasks.push_back(std::move(batch)); }
                task_cv.notify_one();
            }
            { std::lock_guard<std::mutex> lock(mutex); done=true; }
            task_cv.notify_all();
            while(written<issued) collect();
            for(auto& worker:workers) worker.join();
        } catch(...) {
            { std::lock_guard<std::mutex> lock(mutex); done=true; }
            task_cv.notify_all();
            for(auto& worker:workers) if(worker.joinable()) worker.join();
            throw;
        }
        progress.stop();
        if(total.reads==0) throw std::runtime_error("No FASTA records found");
        output.flush();
        quarantine_ids.flush(); direct_repeat_ids.flush();
        if(retained_fasta.is_open()) retained_fasta.flush();
        if(!output || !quarantine_ids || !direct_repeat_ids ||
           (retained_fasta.is_open() && !retained_fasta)) throw std::runtime_error("Failed writing output");
        const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
        std::ofstream summary(o.summary,std::ios::binary);
        if(!summary) throw std::runtime_error("Cannot open summary: "+o.summary);
        summary<<"{\n  \"input\": ";json_string(summary,o.input);
        summary<<",\n  \"adapters\": ";json_string(summary,o.adapters);
        summary<<",\n  \"retained_fasta\": ";json_string(summary,o.retained_fasta);
        summary<<",\n  \"reads\": "<<total.reads<<",\n  \"bases\": "<<total.bases
               <<",\n  \"exclude_adapter\": "<<total.exclude_adapter
               <<",\n  \"exclude_palindrome\": "<<total.exclude_palindrome
               <<",\n  \"quarantine\": "<<total.quarantine
               <<",\n  \"keep_direct_repeat\": "<<total.keep_direct_repeat
               <<",\n  \"keep\": "<<total.keep
               <<",\n  \"adapter_hits\": "<<total.adapter
               <<",\n  \"direct_repeat_hits\": "<<total.direct
               <<",\n  \"inverted_repeat_hits\": "<<total.inverted
               <<",\n  \"threads\": "<<o.threads<<",\n  \"min_repeat_span\": "<<o.min_repeat_span
               <<",\n  \"min_seed_pairs\": "<<o.min_seed_pairs
               <<",\n  \"seed_mod\": "<<o.seed_mod
               <<",\n  \"seconds\": "<<std::fixed<<std::setprecision(3)<<seconds
               <<",\n  \"reads_per_second\": "<<std::fixed<<std::setprecision(1)
               <<(seconds>0?total.reads/seconds:0)<<",\n  \"fasta_quality_and_bam_tags\": \"unavailable\"\n}\n";
        summary.close();
        if(!summary) throw std::runtime_error("Failed writing summary");
        std::cerr<<"[done] "<<total.reads<<" reads; "<<total.bases<<" bases; elapsed="<<std::fixed<<std::setprecision(1)<<seconds<<"s; quarantined="
                 <<(total.exclude_adapter+total.exclude_palindrome+total.quarantine)
                 <<" direct-repeat-candidates="<<total.keep_direct_repeat<<"; "<<std::fixed<<std::setprecision(1)
                 <<(seconds>0?total.reads/seconds:0)<<" reads/s\n";
    } catch(const std::exception& e) {
        std::cerr<<"Error: "<<e.what()<<'\n';
        return 1;
    }
}
