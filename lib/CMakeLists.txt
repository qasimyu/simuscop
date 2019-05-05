cmake_minimum_required(VERSION 2.8)

# Build the matrix library
include_directories(matrix)

# Build the snp library
include_directories(snp)
add_library(snp snp/snp.cpp)

# Build the vcfparser library
include_directories(vcfparser)
add_library(vcfparser vcfparser/vcfparser.cpp)

# Build the split library
include_directories(split)
add_library(split split/split.cpp)

# Build the threadpool library
include_directories(threadpool)
add_library(threadpool threadpool/ThreadPool.cpp)

# Build the seqwriter library
include_directories(seqwriter)
add_library(seqwriter seqwriter/SeqWriter.cpp)

# Build the mydefine library
include_directories(mydefine)
add_library(mydefine mydefine/MyDefine.cpp)
target_link_libraries(mydefine config genome profile threadpool seqwriter)

# Build the config library
include_directories(config)
add_library(config config/Config.cpp)
target_link_libraries(config mydefine split)

# Build the fastahack library
include_directories(fastahack)
add_library(fastahack fastahack/Fasta.cpp)
target_link_libraries(fastahack split)

# Build the genome library
include_directories(genome)
add_library(genome genome/Genome.cpp)
target_link_libraries(genome segment fastahack snp vcfparser split mydefine)

# Build the profile library
include_directories(profile)
add_library(profile profile/Profile.cpp)
target_link_libraries(profile split mydefine psifunc)

# Build the segment library
include_directories(segment)
add_library(segment segment/Segment.cpp)
target_link_libraries(segment fastahack snp mydefine)

# Build the psifunc library
include_directories(psifunc)
add_library(psifunc psifunc/psiFunc.cpp)
target_link_libraries(psifunc)
