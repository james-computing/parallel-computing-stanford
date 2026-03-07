CXX = g++
CXXFLAGS = -std=c++23 -O3 #-ggdb
LDFLAGS = -mavx2 # for SIMD
ISPC = ispc
ISPCFLAGS = -O2 --target=avx
SRC_DIR = src
SRC_CPP_DIR = $(SRC_DIR)/cpp
SRC_ISPC_DIR = $(SRC_DIR)/ispc
OBJ_DIR = obj
OBJ_CPP_DIR = $(OBJ_DIR)/cpp
OBJ_ISPC_DIR = $(OBJ_DIR)/ispc
BIN_DIR = .
SRC_CPP = $(wildcard $(SRC_CPP_DIR)/*.cpp)
SRC_ISPC = $(wildcard $(SRC_ISPC_DIR)/*.ispc)
OBJ_CPP = $(SRC_CPP:$(SRC_CPP_DIR)/%.cpp=$(OBJ_CPP_DIR)/%.o)
OBJ_ISPC = $(SRC_ISPC:$(SRC_ISPC_DIR)/%.ispc=$(OBJ_ISPC_DIR)/%_ispc.o)
INCLUDE_ISPC_DIR = include/ispc
EXE = $(BIN_DIR)/main

.PHONY: all clean

all: $(EXE)

$(EXE): $(INCLUDE_ISPC_DIR) $(OBJ_ISPC) $(OBJ_CPP) $(BIN_DIR)
	$(CXX) -o $(EXE) $(OBJ_ISPC) $(OBJ_CPP) $(LDFLAGS)

$(OBJ_ISPC_DIR)/%_ispc.o: $(SRC_ISPC_DIR)/%.ispc | $(OBJ_ISPC_DIR)
	$(ISPC) $(ISPCFLAGS) $< -o $@ -h $(subst $(OBJ_ISPC_DIR),$(INCLUDE_ISPC_DIR),$(subst .o,.h,$@))

$(OBJ_CPP_DIR)/%.o: $(SRC_CPP_DIR)/%.cpp | $(OBJ_CPP_DIR)
	$(CXX) $(CXXFLAGS) -c $< $(LDFLAGS) -o $@

$(OBJ_CPP_DIR) $(OBJ_ISPC_DIR): $(OBJ_DIR)
	mkdir -p $@

$(OBJ_DIR) $(BIN_DIR) $(INCLUDE_ISPC_DIR):
	mkdir -p $@

clean:
	rm -r $(EXE) $(OBJ_DIR) $(INCLUDE_ISPC_DIR)

-include: $(OBJ:.o=.d)