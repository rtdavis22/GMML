// Author: Robert Davis

#ifndef GMML_INTERNAL_STUBS_FILE_H_
#define GMML_INTERNAL_STUBS_FILE_H_

#include <iosfwd>
#include <stdexcept>
#include <string>

namespace gmml {

// This is similar to the File class in the Java API, but way less robust.
class File {
  public:
    explicit File(const std::string& pathname) : pathname_(pathname) {}

    File(const File& dir, const std::string& pathname)
            : pathname_(dir.pathname() + pathname) {}

    bool exists() const;

    void set_pathname(const std::string& pathname) { pathname_ = pathname; }

    std::string pathname() const { return pathname_; }

    const char *c_str() const { return pathname_.c_str(); }

  private:
    std::string pathname_;
};

class FileNotFoundException : public std::invalid_argument {
  public:
    explicit FileNotFoundException(const std::string& file_name)
            : std::invalid_argument("File not found: " + file_name),
              file_name_(file_name) {}

    ~FileNotFoundException() throw() {}

    std::string file_name() const { return file_name_; }

  private:
    std::string file_name_;
};

class Readable {
  public:
    // This should be removed.
    void read(const std::string& file_name);

    // Throws FileNotFoundException.
    void read(const File& file);

  private:
    virtual void read(std::istream&) = 0;
};

class Writeable {
  public:
    void print() const;

    // This should be removed.
    void print(const std::string& file_name) const;

    void print(const File& file) const;

  private:
    virtual void write(std::ostream&) const = 0;
};

}  // namespace gmml

#endif  // GMML_INTERNAL_STUBS_FILE_H_
