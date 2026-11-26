########################################################################
#
# This file is part of OpenPMX (https://github.com/deleveld/openpmx).
# Copyright (c) 2024 Douglas Eleveld.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

# OpenPMX Makefile (replaces the install script)
#
# usage: make [PARALLEL=pthread|openmp|singlethread] [SERVER=notserver|server]
#        make clean
#        make -C /path/to/openpmx ...
#
# Changing any setting (or the compiler) rebuilds everything that depends on
# it. PLATFORM (macos, linux, wsl, msys2) is auto-detected and can be
# overridden.

OPENPMXPATH := $(CURDIR)

PARALLEL ?= pthread
SERVER   ?= notserver

# compiler, honour CC from the environment or command line
ifeq ($(origin CC),default)
CC := gcc
endif

# find out about the platform
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
DETECTED := macos
else ifeq ($(UNAME_S),Linux)
DETECTED := linux
else ifneq (,$(filter MSYS% MINGW% CYGWIN%,$(UNAME_S)))
DETECTED := msys2
else
DETECTED := unknown
endif
PLATFORM ?= $(DETECTED)
ifeq (,$(filter macos linux wsl msys2,$(PLATFORM)))
$(error platform $(PLATFORM) not recognized)
endif

# Windows executables need .exe
EXEEXT :=
ifeq ($(PLATFORM),msys2)
EXEEXT := .exe
endif

# general build flags
LIBNAME   := openpmx
OPX_CFLAGS    := -g -O2 -W -Wall -Wextra -march=native -mtune=native
OPX_CLIBFLAGS := -Wfloat-conversion -Wformat -Werror=format-security -Wno-override-init -Wwrite-strings -Isrc
OPX_LFLAGS    :=
OPX_LDLIBS    := -lgsl -lgslcblas -lm
LAUNCHER_CFLAGS := -g -O2 -Wall

# adjust flags for different platforms
ifeq ($(PLATFORM),macos)
OPX_CLIBFLAGS += -fPIC
ifeq ($(shell uname -m),arm64)
OPX_CFLAGS += -mcpu=apple-m1
else
OPX_CFLAGS += -march=native
endif
OPX_LFLAGS += -flto
else ifneq (,$(filter linux wsl,$(PLATFORM)))
OPX_CLIBFLAGS += -fPIC
OPX_LFLAGS += -flto
else ifeq ($(PLATFORM),msys2)
OPX_CFLAGS += -I/mingw64/include
OPX_LFLAGS += -L"$(OPENPMXPATH)/lib" -L/mingw64/lib
OPX_LFLAGS += -flto
endif

# handle different flags for parallel targets
ifeq ($(PARALLEL),openmp)
OPX_CFLAGS += -fopenmp
OPX_LFLAGS += -fopenmp
OPX_LDLIBS += -lgomp
PARALLEL_DEFINE := OPENPMX_PARALLEL_OPENMP
else ifeq ($(PARALLEL),pthread)
OPX_CFLAGS += -pthread
OPX_LFLAGS += -pthread
PARALLEL_DEFINE := OPENPMX_PARALLEL_PTHREADS
else ifeq ($(PARALLEL),singlethread)
PARALLEL_DEFINE := OPENPMX_PARALLEL_SINGLETHREAD
else
$(error PARALLEL "$(PARALLEL)" not recognized (pthread, openmp, singlethread))
endif

# for server the queue is in the install dir so only a single instance
# can run per install.
ifeq ($(SERVER),server)
SERVER_DEFINE := OPENPMX_SERVER
else ifeq ($(SERVER),notserver)
SERVER_DEFINE := OPENPMX_NOTSERVER
else
$(error SERVER "$(SERVER)" not recognized (server, notserver))
endif

# quote for the shell, and escape for embedding in a C string literal
shq  = '$(subst ','\'',$1)'
cstr = "$(subst ",\",$(subst \,\\,$1))"

# source code, object files are flat in build/
SOURCEFILES := pmxstate.c options.c ievaluate.c predict.c stage1.c \
	checkout.c popmodel.c encode.c nonzero.c omegainfo.c \
	linalg.c print.c idata.c table.c scatter.c estimate.c \
	simulate.c reload.c server_queue.c \
	profile.c covariance.c vpc.c npde.c \
	dataconfig/dataconfig.c \
	advan/advan.c advan/pred.c \
	advan/onecomp.c advan/onecomp_depot.c \
	advan/twocomp.c advan/threecomp.c \
	advan/diffeqn_test.c advan/diffeqn_libgsl.c \
	advan/eigen.c advan/tci.c \
	models/propofol_eleveld.c \
	models/propofol_schnider.c \
	models/remifentanil_eleveld.c \
	models/remifentanil_minto.c \
	models/sufentanil_gepts.c \
	models/remimazolam_eleveld.c \
	utils/vector.c utils/various.c \
	utils/getdelim.c utils/errctx.c utils/filecontent.c \
	bobyqa/bobyqa.c
OBJECTS := $(addprefix build/,$(addsuffix .o,$(basename $(notdir $(SOURCEFILES)))))
vpath %.c $(sort $(addprefix src/,$(dir $(SOURCEFILES))))
vpath %.c src

TRAN     := bin/openpmxtran$(EXEEXT)
LIB      := lib/lib$(LIBNAME).a
LAUNCHER := openpmx$(EXEEXT)
BUILDFLAGS := src/buildflags.h

# everything that influences the build, rewritten only when it changes
CONFIG := $(PLATFORM) $(PARALLEL) $(SERVER) $(CC) $(OPX_CFLAGS) $(OPX_CLIBFLAGS) $(OPX_LFLAGS) $(OPX_LDLIBS) $(LAUNCHER_CFLAGS) $(OPENPMXPATH)

.PHONY: all clean help FORCE
all: $(TRAN) $(LIB) $(LAUNCHER)
	@echo "openpmx install success ($(PARALLEL) $(SERVER) $(PLATFORM))"

help:
	@echo "This file installs OpenPMX (https://github.com/deleveld/openpmx)."
	@echo "usage: make [PARALLEL=pthread|openmp|singlethread] [SERVER=notserver|server] | make clean"

build bin lib:
	mkdir -p $@

build/config: FORCE | build
	@printf '%s\n' $(call shq,$(CONFIG)) > $@.new; \
	if cmp -s $@.new $@; then rm -f $@.new; else mv $@.new $@; fi

$(BUILDFLAGS): build/config
	@echo "generate $@"
	@{ printf '#define OPENPMX_INSTALL_PLATFORM_%s\n' $(call shq,$(PLATFORM)); \
	   printf '#define OPENPMX_INSTALL_PLATFORM_NAME "%s"\n' $(call shq,$(PLATFORM)); \
	   printf '#define OPENPMX_INSTALL_PATH %s\n' $(call shq,$(call cstr,$(OPENPMXPATH))); \
	   printf '#define %s\n' $(PARALLEL_DEFINE); \
	   printf '#define %s\n' $(SERVER_DEFINE); \
	} > $@

# compile openpmxtran
$(TRAN): src/openpmxtran.c $(BUILDFLAGS) build/config | bin
	@echo "compile openpmxtran"
	$(CC) $(OPX_CFLAGS) src/openpmxtran.c -o $@ $(OPX_LFLAGS) -lm

# compile source code to object files
$(OBJECTS): $(BUILDFLAGS) build/config
build/%.o: %.c | build
	@echo "compile $(patsubst src/%,%,$<)"
	$(CC) $(OPX_CFLAGS) $(OPX_CLIBFLAGS) -Iinclude -MMD -MP -c $< -o $@

# build library from object files (start fresh so no stale members remain)
$(LIB): $(OBJECTS) | lib
	@echo "build library"
	rm -f $@
	ar rcs $@ $^

# compile the C launcher
$(LAUNCHER): src/openpmx_launcher.c src/utils/c22.h build/config
	@echo "compile openpmx launcher"
	$(CC) $(LAUNCHER_CFLAGS) \
		$(call shq,-DOPENPMX_CC=$(call cstr,$(CC))) \
		$(call shq,-DOPENPMX_CFLAGS=$(call cstr,$(OPX_CFLAGS))) \
		$(call shq,-DOPENPMX_LFLAGS=$(call cstr,$(OPX_LFLAGS))) \
		$(call shq,-DOPENPMX_LDLIBS=$(call cstr,$(OPX_LDLIBS))) \
		$(call shq,-DOPENPMX_PATH=$(call cstr,$(OPENPMXPATH))) \
		src/openpmx_launcher.c -o $@

clean:
	rm -vf build/*.o build/*.d build/config build/config.new
	rm -vf bin/openpmxtran bin/openpmxtran.exe
	rm -vf lib/*.a
	rm -vf openpmx openpmx.exe
	rm -vf $(BUILDFLAGS)
	@echo "openpmx clean success"

-include $(OBJECTS:.o=.d)
