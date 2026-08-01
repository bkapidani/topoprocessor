#include "mapped_file.h"

#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>

namespace {

void require(bool condition)
{
    if (!condition) {
        throw std::runtime_error("mapped_file test assertion failed");
    }
}

class temporary_files {
public:
    explicit temporary_files(const std::string& directory)
        : text_(directory + "/mapped_file_test.txt"),
          empty_(directory + "/mapped_file_empty.txt")
    {}

    ~temporary_files()
    {
        std::remove(text_.c_str());
        std::remove(empty_.c_str());
    }

    const std::string& text() const { return text_; }
    const std::string& empty() const { return empty_; }

private:
    std::string text_;
    std::string empty_;
};

void write_file(const std::string& filename, const std::string& contents)
{
    std::ofstream output(filename, std::ios::binary);
    require(static_cast<bool>(output));
    output.write(contents.data(), static_cast<std::streamsize>(contents.size()));
    require(static_cast<bool>(output));
}

} // namespace

int main(int argc, char** argv)
{
    require(argc == 2);
    temporary_files files(argv[1]);
    write_file(files.text(), "first\nsecond\nthird");
    write_file(files.empty(), "");

    mapped_file input(files.text());
    require(input.is_open());
    require(!input.end());
    require(input.mem() != nullptr);
    require(std::string(input.mem(), 5) == "first");
    require(input.mem()[18] == '\0');
    require(input.get_line() == "first");
    require(input.get_line() == "second");
    require(input.get_line() == "third");
    require(input.end());

    bool exhausted = false;
    try {
        input.get_line();
    } catch (const std::out_of_range&) {
        exhausted = true;
    }
    require(exhausted);

    input.open(files.empty());
    require(input.is_open());
    require(input.end());
    require(input.mem() != nullptr);
    require(*input.mem() == '\0');

    input.open(files.text() + ".missing");
    require(!input.is_open());
    require(input.end());
    require(input.mem() == nullptr);

    input.open(files.text());
    require(input.is_open());
    input.close();
    require(!input.is_open());
    require(input.end());
    require(input.mem() == nullptr);

    return 0;
}
