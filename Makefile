CXX ?= g++
CXXFLAGS ?= -O3 -Wall -Wextra
PREFIX ?= $(HOME)/.local

# GCC 8 ships C++17 filesystem in a separate library. Probe the selected
# compiler and flags, rather than assuming a compiler version or vendor.
# Override with FS_LIBS=-lstdc++fs (or FS_LIBS=) if cross-compiling.
FS_LIBS = $(shell printf '%s\n' '\#include <filesystem>' 'int main() { return std::filesystem::exists("."); }' | $(CXX) $(CPPFLAGS) $(CXXFLAGS) -std=c++17 -pthread -x c++ - -o /dev/null $(LDFLAGS) $(LDLIBS) >/dev/null 2>&1 || printf '%s' '-lstdc++fs')

.PHONY: all test install clean

all: bin/hificheck

bin/hificheck: hificheck.cpp Makefile
	mkdir -p bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -std=c++17 -pthread -o $@ $< $(LDFLAGS) $(LDLIBS) $(FS_LIBS)

test: bin/hificheck
	python3 -m unittest discover -s tests -p 'test_hificheck.py' -v

install: bin/hificheck
	install -d $(DESTDIR)$(PREFIX)/bin
	install -m 755 bin/hificheck $(DESTDIR)$(PREFIX)/bin/hificheck

clean:
	rm -f bin/hificheck
