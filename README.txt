/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * Copyright (C) 2024-2025 Ramanakumar Sankar                      *
 * Copyright (C) 1998-2023 Timothy E. Dowling                      *
 *                                                                 *
 * This program is free software; you can redistribute it and/or   *
 * modify it under the terms of the GNU General Public License     *
 * as published by the Free Software Foundation; either version 2  *
 * of the License, or (at your option) any later version.          *
 * A copy of this License is in the file:                          *
 *   $EPIC_PATH/License.txt                                        *
 *                                                                 *
 * This program is distributed in the hope that it will be useful, *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of  *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.            *
 *                                                                 *
 * You should have received a copy of the GNU General Public       *
 * License along with this program; if not, write to the Free      *
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,     *
 * Boston, MA 02110-1301, USA.                                     *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

 Explicit Planetary Isentropic-Coordinate 
      (EPIC) Atmospheric Model
            Version 5.30

  Required: 
    Unix operating system
    C compiler
    CMake, for managing source-code compilation and linking
    netCDF, for reading and writing self-describing .nc files

  Optional: 
    MPI (Message Passing Interface), for running on parallel computers
    
    A host of software packages will plot netcdf (.nc) files.
    A good one to start with is Panoply, which is free from NASA GISS:
      http://www.giss.nasa.gov/tools/panoply/

    Sample tools for Matlab, IDL, Mathematica, etc. are included in
    $EPIC_PATH/tools, and additions are encouraged. 

 Topics covered:
   A. INSTALLING
   B. TROUBLESHOOTING
   C. INITIALIZING
   D. CHANGING
   E. RUNNING FROM COMMAND-LINE
    
=========================================================================
=========================================================================

A. INSTALLING 

  1.  To decompress and de-archive the file epic.tar.gz:

      a) If a directory named epic already exists, delete it (or move it):
           rm -r epic
      b) Type 
           tar zxpf epic.tar.gz

      This will produce a directory tree similar to:

%cd $EPIC_PATH
%ls -l -a

epic/             Top EPIC-model directory
  bin/            Executable files
  CMakeLists.txt  The script that cmake runs to compile the source code
  data/           Data files (planetary data, chemical data, etc.)
  help/           Help files for executables
  include/        Header (.h) files
  License.txt     Software license
  notes/          Brief remarks on version history, etc.
  README.txt      This file.
  src/
    attic/        A place to store valuable old code.
    chemistry/    Subroutines related to atmospheric chemistry
    clouds/       Cloud microphysics code
    core/         Dynamical core code
    mpi/          Parallel-processor specific code
    rt/           Radiative transfer code
    single/       Single-processor specific code
    subgrid/      Turbulence and other subgrid code
  tools/          Analysis tools that work with external software
  util/           Unix utility programs such as shell scripts

  2. To archive and compress the model:

      a) Clear any object code by typing
           cd $EPIC_PATH/build
           make clean

      b) Move up to the directory containing epic and type
           touch epic
           tar cf epic.tar epic
           gzip epic.tar

     Note: usage of the environment variable $EPIC_PATH is enabled after completing step 3.


  3. In your shell resource file (such as .cshrc), add the
     following environment-variable lines, and then edit them for your environment.

     For a .cshrc or .tcshrc file:

setenv CC        clang
setenv EPIC_PATH ~/epic

      For a .bashrc file:

export CC=clang
export EPIC_PATH=~/epic

      EPIC_PATH is the unix path where the EPIC model is kept.
      CMAKE_BUILD_TYPE can be set to Release, Debug, or RelWithDebInfo

      To set the environment variables in the .cshrc resource file for the current shell, type
        source ~/.cshrc
                 
  4.  To compile the model:

      a) To compile the model the first time:
           cd $EPIC_PATH
           mkdir build
           cd build
           cmake -DCMAKE_BUILD_TYPE=Release ..
           make install

         These steps do not need to be repeated unless the platform's netCDF
         or MPI software changes, or a different CMAKE_BUILD_TYPE is desired
         (the options are Release, Debug, and RelWithDebInfo). 

      b) In general (after the first time), to compile, link, and install:
           cd $EPIC_PATH/build
           make install

      c) Prior to archiving, pare down and keep just the source code:
           cd $EPIC_PATH/build
           make clean
           cd ..
           rm -r -f build

      The executables are installed in $EPIC_PATH/bin.

=========================================================================
=========================================================================

B. INITIALIZING

   The information needed to start the EPIC model is contained in
   a single epic.nc file. The program "initial" generates this file.

=========================================================================
=========================================================================

C. CHANGING

   To make select changes to parameters of in an existing epic.nc file,
   use the program "change," which reads an epic.nc file, prompts for
   changes, and then writes a new file. You may wish to customize the
   prompts in $EPIC_PATH/src/shared/epic_change.c. 

   To learn about other uses of the change program, include the "-h" or
   "-help" command-line flag.

=========================================================================
=========================================================================

D. RUNNING FROM COMMAND-LINE

   Typing the name of the epic program with the command-line flag
   -h or -help prints the options on the screen.  For example, type:

%epic -h

   to get a list of command-line flags and their meanings.
      
   To integrate the file epic.nc forward 10000 timesteps as a batch job,
   saving twice and backing up every 1000, with error messages to be 
   written to the file ``epic.log'', type:

%epic -itrun 10000 -itsave 5000 -itback 1000 epic.nc >& epic.log &

   To write an extract.nc file every 100 steps, first choose the variables
   to extract using initial or change, and then run the model with:

%epic -itrun 10000 -itsave 5000 -itback 1000 -itextract 100 epic.nc >& epic.log &

   To append to an existing extract.nc file called extract2.nc, type

%epic -itrun 10000 -itsave 5000 -itback 1000 -itextract 100 -append extract2.nc epic.nc >& epic.log &

=========================================================================
=========================================================================

