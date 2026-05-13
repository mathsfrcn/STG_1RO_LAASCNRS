
#SYSTEM     = x86-64_osx
SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

CC = g++ -O0
CFLAGS=-g -W -m64 -O -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD -std=c++17
#CFLAGS=-W -Wall -m64 -O -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD 

#CPLEX_DIR=/Users/idir/Applications/IBM/ILOG/CPLEX_Studio1263
CPLEX_DIR=/net/cubitus/projects/Partage_ROC/Logiciels/CPLEX_Studio2212
INCLUDES=-I$(CPLEX_DIR)/cplex/include -I$(CPLEX_DIR)/concert/include
LIBDIRS=-L$(CPLEX_DIR)/cplex/lib/x86-64_linux/static_pic/ -L$(CPLEX_DIR)/concert/lib/x86-64_linux/static_pic/
LIBS=-lconcert -lilocplex -lcplex -lm -lpthread -ldl 
#Valable pour le 1271
#LIBS=-lilocplex -lcplex -lconcert -lm -lpthread 
SRCS=$(wildcard *.cpp)
OBJS=$(SRCS:.cpp=.o)
EXEC=benders_kc

all: $(EXEC) clean

$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) -o $(EXEC) $(OBJS) $(LIBDIRS) $(LIBS) 

%.o: %.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	rm -f *.o