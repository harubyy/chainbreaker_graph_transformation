CXX      := $(CXX)
CXXFLAGS := -O3 -std=c++17 -fopenmp
LDFLAGS  := -lufget -lmatio 
BUILD    := ./build
OBJ_DIR  := $(BUILD)/objects
APP_DIR  := $(BUILD)
#TARGET   := rewriting
INCLUDE  := -Iinclude/
SRC      :=                      \
   $(wildcard src/*.cpp)         \

OBJECTS := $(SRC:%.cpp=$(OBJ_DIR)/%.o)

#all: build $(APP_DIR)/$(TARGET)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

$(APP_DIR)/release: $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(OBJECTS) $(LDFLAGS) -o $(APP_DIR)/release

$(APP_DIR)/rewrite: $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(OBJECTS) $(LDFLAGS) -o $(APP_DIR)/rewrite

$(APP_DIR)/rewrite_dump: $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(OBJECTS) $(LDFLAGS) -o $(APP_DIR)/rewrite_dump

$(APP_DIR)/debug: $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(OBJECTS) $(LDFLAGS) -o $(APP_DIR)/debug

.PHONY: build clean debug release rewrite rewrite_dump

build:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)

debug:		CXXFLAGS += -DDEBUG -DREWRITE_ENABLED -g
debug:		build $(APP_DIR)/debug

release:	build $(APP_DIR)/release
	
rewrite:	CXXFLAGS += -DREWRITE_ENABLED
rewrite:	build $(APP_DIR)/rewrite
	
#rewrite_dump:	CXXFLAGS += -DREWRITE_ENABLED -DDUMP_CSR
#rewrite_dump:	build $(APP_DIR)/rewrite_dump

clean:
	-@rm -rvf $(OBJ_DIR)/*
#	-@rm -rvf $(APP_DIR)/*
