CXX ?= g++
CXXFLAGS ?= -O3 -std=c++17 -pthread -Wall -Wextra
PREFIX ?= $(HOME)/.local

.PHONY: all test install clean

all: bin/hificheck

bin/hificheck: hificheck.cpp
	mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ $<

test: bin/hificheck
	python3 -m unittest discover -s tests -p 'test_hificheck.py' -v

install: bin/hificheck
	install -d $(DESTDIR)$(PREFIX)/bin
	install -m 755 bin/hificheck $(DESTDIR)$(PREFIX)/bin/hificheck

clean:
	rm -f bin/hificheck
