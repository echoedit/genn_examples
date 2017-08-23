ifdef MPI_ENABLE
EXECUTABLE      := simulator_${OMPI_COMM_WORLD_RANK}
else
EXECUTABLE      := simulator
endif
SOURCES         := simulator.cu
include $(GENN_PATH)/userproject/include/makefile_common_gnu.mk
