
CXX=g++
CXXFLAGS = --std=c++17 -O3 -march=native
CXXFLAGS += -I/home/tomgrad/include

qvoter : main.o
	$(CXX) -lpthread  -o $@ $^

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -rf *.o qvoter

