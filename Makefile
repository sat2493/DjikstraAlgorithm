CXXFLAGS = -std=c++11 -Wall -Werror

all: shortest_path

shortest_path: shortest_path.cc
	$(CXX) -o $@ $^

execs = shortest_path

clean:
	rm -f $(execs)

