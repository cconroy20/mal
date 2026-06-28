                COWAN'S CODE MODIFICATION FOR IBM PC
                            by A. KRAMIDA
                (Adapted also for IBM RISC/6000 in NIST, USA,
                VAC-VMS in Bohum, Germany, and other platforms)

        National Institute of Standards and Technology,
                    Gaithersburg, MD 20899, USA
                e-mail: alexander.kramida @ nist.gov
                       tel. +1 (301) 975-8074

                       Distributed from
     https://datapub.nist.gov/od/id/6CF509047B474AC9E05324570681DE731930

             Also distributed by the Institute for Spectroscopy,
                     Russian Academy of Sciences,
                TROITSK, MOSCOW REGION 142092 RUSSIA
                http://das101.isan.troitsk.ru/COWAN

        Current version release date: November 1, 2021

CONTENTS
        The package consists of a set of computer programs for calculation
of energy levels, radiative transition wavelengths and probabilities,
electron impact excitation and photoionization cross sections etc. The
programs were originally written by Robert D. Cowan in 1962-1991
(see R.D. Cowan, The theory  of atomic structure and spectra, University
of California Press, Berkeley - Los Angeles - London, 1981).
Several people contributed to their texts. Complete instructions on
how to use these programs are available in the documentation (README)
files coming along with this package. You can contact A. Kramida by
e-mail or telephone given above if you have any questions on the programs
operation. The latest versions of the original R. Cowan's codes are
available at

        http://www.tcd.ie/Physics/People/Cormac.McGuinness/Cowan/

        In addition to the four main programs, RCN, RCN2, RCG and RCE
(there are both executable files and Fortran source files), the package
contains many supplementary files. Some of them are used by the main
programs, and others are helpful utilities. If you want the package to
operate properly, you should follow the installation procedure described
below. The complete list of all files in the package is given in the
addendum below.
        The package contains executable files that are configured to be
used with tasks that fit within definite array dimensions. If you need to
use it for larger tasks, you should re-compile the programs and (maybe)
re-build some data files used by RCG (instructions are in the RCG
manual). For PCs, we recommend to use the gfortran compiler within the
MinGW distribution, which is freely available for non-commercial use and 
compatible with the Windows operating system. The codes can also be 
compiled on various other computer platforms. For details, see section
COMPILATION IN OTHER OPERATING SYSTEMS below.

CURRENT VERSIONS
RCN  ver. 36k, March 2019,
RCN2 ver. 2k,  March 2019,
RCG  ver. 11k, April 2021,
RCE  ver. 20k, September 2021.

INSTALLATION (Only for Windows-based PCs; mainframe and workstation users
can skip this section)

        The installation procedure is designed for installation from a
compressed (ZIP) file downloaded in its entirety from the Web.

        In order to install the package on your PC, you should make the
following steps:
        1) Extract all files from the COWAN.ZIP archive PRESERVING THE
SUB-DIRECTORY STRUCTURE into the directory that you choose to host the
package (typically, c:\COWAN; however, you can use whatever destination
is suitable for you).
        2) Modify the system's environment. Namely, the path to the Cowan
executable files has to be appended to the system PATH variable. If you
are installing the package in a network location, you may need to
assign a drive letter to it. Windows command line interpreter (cmd.exe)
does not support network paths (starting with '\\'). Ask your system
administrator for assistance.
        Windows XP: Right-click on the 'My Computer' icon on the desktop
(if you do not have this icon on the desktop, right-click on an empty
space on the desktop, select 'Properties' in the pop-up menu, select tab
'Desktop', click 'Customize Desktop' and check the box 'My Computer').
In the pop-up menu, select 'Properties'. Select tab 'Advanced' and click
the button 'Environment variables'. There are two panes: User-defined
variables and System variables. You can edit the PATH in any of them.
The changes you make in the System Variables pane are effective for all
users, while the changes in the User Variables pane are effective only
for the current user. Do the change in the System Variables pane if you
have sufficient privileges. Otherwise, do them in the User Variables
pane. Select the PATH variable in the list and click Edit. Append the
path string with the path to the COWAN\CODE directory, e.g.,
...;c:\COWAN\CODE
where ... denotes the previous contents of the PATH string. Click Ok to
close the editing dialog, Ok to close the Environment Variables box, and
Ok to close the System Properties box.
        Windows 7 and other versions of Windows later than XP: The
'Environment Variables' settings are accessible through the Control Panel.
Consult Windows Help for more details on how to change the PATH variable.
        3) Open the command prompt window. On Windows XP and Windows 7, there
are two ways of doing it: 1) Start/Command Prompt; 2) Run/cmd (type 'cmd'
in the prompt box and hit <Enter>).
        4) Adjust the width of the Command Prompt window to at least 120
characters (right-click on the window title, select Properties, change
Window Size Width, make any other adjustments you want, e.g. font size,
colors, etc., click Ok, select 'Save properties for future windows with
same title', click Ok.
        5) Type PATH and hit 'Enter'. If you do not see the path to the
COWAN\CODE directory in the PATH string, close the Command Prompt window
(either click on the cross in the right upper corner or type 'Exit' and
hit 'Enter'), and return to Step 2.
        6) Change current directory to COWAN\WORK, e.g.
cd c:\COWAN\WORK <Enter>.
Execute the make_cfp.bat batch file:
        make_cfp <Enter>
This builds the binary CFP 'decks' (files FOR072, FOR073, FOR074 in
COWAN\CODE). It should take about half a minute, and you should receive
a message 'Done. New CFP binary files for072, ... are created in
...COWAN\CODE\ directory'. If you do not see this message, or any of the
files FOR072, FOR073, FOR074 has a zero size, RCG will not work! Contact
A. Kramida to resolve the problem.
        7) Execute the test task defined by the sample input files IN36 and
IN2 in the COWAN\WORK directory:
        RCN  <Enter>
        RCN2 <Enter>
        RCG  <Enter>
        RCE  <Enter>
Each of these commands should report a 'NORMAL EXIT'. If any of them
fails, installation was not successful. Contact A. Kramida to resolve
the problems.

        After successful installation, the directory COWAN should contain
this README_Kramida.txt file, Cowan's original manuals Readme.cowan.htm,
RCN_DOC.txt, RCG_DOC.txt, and RCE_DOC.txt and the following
sub-directories:
        CODE - this directory should contain all executable codes and data
        files needed for them.
        FOR  - this directory will contain the Fortran source texts of the
        main programs.
        WORK - this directory will contain sample input files IN36 and
        IN2. Starting with these files, you can run all the
        programs subsequently to make a test that the installation
        was completed properly. The files IN36.ai and IN2.ai are
        sample input files for autoionization calculations.

In addition to the main programs RCN, RCN2, RCG, and RCE, the package
includes a number of utility codes described below. Most of them are
written in Perl and require Perl to be installed. I recommend
installing the Strawberry Perl, which can easily be found on the
Internet. It is freely available under the GNU public license. Some of
the old utility codes written and compiled using Borland Pascal will
not run under Windows 7. However, they still can be used with the
freely available DosBox software.


RUNNING THE PROGRAMS
        Programs are normally run in the following order:
RCN, RCN2, RCG, RCE.
        For PCs, these names are names of batch files that must be entered
at the command prompt from the directory containing the input files.
Each program produces the input files for the next one. After the first
run, you may have multiple runs of RCG and/or RCE.
        If Perl is installed, the Perl utilities can be used. The most
commonly used ones are add_eav, update11, printout, and conv_out. See 
their descriptions below.

ADDITIONS TO THE MANUAL.

        1) The input/output file names are given in the table below.
------------------------------------------------------------------------
Program Input file(s)    Output file(s)  Contents of the output file(s)
------------------------------------------------------------------------
RCN    IN36    (ASCII)  OUT36   (ASCII) Diagnostics
                         TAPE2N (binary) Input for RCN2 (wavefunctions)

RCN2   IN2     (ASCII)  OUT2    (ASCII) Diagnostics
        TAPE2N (binary)  ING11   (ASCII) INPUT for RCG  (can be edited)

RCG    ING11   (ASCII)  OUTG11  (ASCII) Diagnostics
                         OUTGINE (ASCII) Input for RCE (can be edited)
                         TAPE2E (binary) Input for RCE

RCE    OUTGINE (ASCII)  OUTE    (ASCII) Diagnostics
        TAPE2E (binary)  LEVELS1 (ASCII)   -"-
        [RCEINP] (ASCII)  LEVELS2 (ASCII)   -"-
                         LEVELS3 (ASCII)   -"-
                         PARVALS (ASCII)   -"-
                         RCEOUT  (ASCII) Can be renamed to and used as
                                         RCEINP file
                        [RCEINP] (ASCII) Vector input file for RCE
                     [RCEINP.HF] (ASCII) Copy of RCEINP created first
                                         time you run RCE.
------------------------------------------------------------------------
        File RCEINP can be created by RCE or can be used by RCE as an
input file, depending on a flag INVECT in the OUTGINE file (see later).

        2) Compared to the original Cowan's version, the RCG program
has two more input parameters that must be specified in the
"first control card" (usually, the first line in the file ING11).
They both are in the format INTEGER*1 and are specified in
columns 77 and 78 (instead of disabled parameter ICTC, see the manual
for RCG). Their default values are 0,0. This corresponds to the usual
printout format of RCG, as designed by R. Cowan.  If the second
parameter IQS (column 78 in ING11) is set to non-zero value n<10,
then the quantum numbers that refer to n-th open subshell will appear
in parent term names in the output (this is significant in the case
when there are 3 or more open shells, when it was not obvious to
understand to what shell the numbers in the original Cowan's output
refer). The first of the two parameters, INWOL, does not make any
difference.
        If IQS equals zero and INWOL (column 77 in ING11) equals 1, then
the alternate way of choice of the parent subshell is used. For
questions, refer to the author of these two parameters, Yu. Ralchenko
(present e-mail address yuri.ralchenko @ nist.gov).
        3) Parameter DMIN (columns 61-65 in the RCG input file ING11)
originally could be non-negative floating-point number less than 1.0.
Its meaning was the minimum line strength that will be retained in the
RCG output transitions list. Lacking in the manual, there is also
another way to use this parameter: if it is greater than 1.0, then only
the lines with gA>=10^DMIN sec-1 will be included in the transition
list. This option now works well for any way of line sorting
(originally it did not work at all for wavelength-sorted list).
        4) Designations of terms has slightly changed. Sometimes it
happened (in original version of RCG) that different terms had equal
names, even if the additional (e.g., seniority) quantum numbers are
used. To eliminate this ambiguity in term notation, I introduced for
such terms additional letter ('a', 'b'...) making the final term names
different. For example, the term names (3P) 3Da and (3P) 3Db mean
different 3D terms: their parent terms have the same name but different
genealogy. For such cases, the user may have to manually trace the
history of the summation of quantum numbers that produced these terms,
using some additional printout in OUTG11. This task is automated to
some extent by Perl utility codes printout.bat and conv_out.bat (see
below).
        5) Notation of the configuration-interaction parameters has been
changed. The first two symbols are the sequential numbers of the two
interacting configurations (the same as it was in Cowan's original RCG),
but now the letters 'A','B', ..., 'Z', 'a', ..., 'z' are used for
numbers 10, 11, ... 35, 36, ..., 61, instead of an asterisk ('*')
appeared as FORTRAN replacement for output that does not fit the format
specified (I1). In the latest version of Cowan's RCG, the format of
the CI parameters was extended from 8 to 10 characters, eliminating the
need for such encoding. However, I did not implement this extension
in the present version of the package.
        6) The alternate-coupling scheme input files for RCE (other than
LS/JJ) are created by RCG if the following conditions are fulfilled:
        1. Some non-zero NLSMAX (e.g. 999) is specified in columns 23-25
        in the first line of ING11 (default is blank!).
        2. Only one zero is present in columns 33-37 of the first line
        of ING11 (KCPLD(3..7): no output if non-zero).
        3. If coupling schemes No. 5, 6 or 7 are needed, the number of
        open subshells NOSUBC must be not less than 3 (add more
        open inner shell configurations in IN36 if it is not so!).
        7) Two new optional input parameters have been introduced:
E0MAX(1) and E0MAX(2), which restrict the list of lines included by RCG
in OUTG11 so that the energy of the first-parity levels does not exceed
E0MAX(1) (if it is non-blank), and the energy of the second-parity
levels does not exceed E0MAX(2) (if it is non-blank). These parameters
can be specified in the optional 'rescale card', e.g., on the
additional line included at the top of ING11 with zero in position 5
and floating-point numbers E0MAX(1) and E0MAX(2) in positions 66-75 and
76-85, respectively. The values are given in units of kK (10^3 cm^-1).
        8) The parameter ICRIT2 is disabled in RCE20k. That means, its value
has no effect on the execution of the program: program exits after
completing the fitting either with linked or with unlinked parameters,
and it never starts another run with next value of LSQM.
        9) 'Suspended' execution is not supported in RCE20k: you cannot
save intermediate results after n-th iteration and then start again
from that point. However, the use of the new RCEOUT/RCEINP input files
(see below) permits doing almost the same: the only extra time is spent
on performing an extra initial iteration when making a new run with data
from RCEOUT copied into RCEINP file. The I216 parameter is disabled in
RCE, and its position in OUTGINE file is used by the IW6 parameter
(directing some outpu2:23 PM 8/20/2012t to the screen if IW6<0).
10) In some rare cases when all energy levels are flagged as
"experimental", and all the Slater parameters are set free (varied),
iterations will not be performed by RCE unless the number of iterations
is set to be different from RCG's default value of '5'.
11) Program RCE now permits the use of an alternate input file
format. After the first run, RCE produces new input file called RCEINP.
Apart from it, RCE always creates the file RCEOUT having the same
structure as RCEINP. It makes much easier changing the input data, e.g.
fixing or unfixing energy levels and Slater parameters. File RCEOUT can
be modified with any ASCII editor and renamed as RCEINP. Then it becomes
the new input file for subsequent iterations. A user can ignore these
new features and continue using the old OUTGINE file if so desired.
However, if one wants to use RCEINP file instead of OUTGINE, he/she must
specify the integer flag "1" in the 25th column of the first line in
OUTGINE (parameter INVECT: "use vector-input file"). Then all subsequent
lines in OUTGINE (except the first control line of the next parity
section), will be ignored, and the corresponding data will be read from
RCEINP file.
A new algorithm of sorting the eigenvalues is implemented in RCE. It
has much better performance than Cowan's original options in the sense
of keeping the correspondence between level flags during the iterations,
but it gives no guarantee of keeping the true level flags in all
possible cases. In some cases (probably in situations where the coupling
scheme is very far from pure LS), the user may have a need to increase
the number of eigenvector components used to recognize the energy
levels. The default number is 5. If you see, after some iteration, that
the level flags become incorrect, you may try to increase this value by
specifying a greater parameter MAXC in the source code and recompiling
the program. The same can be done if RCE abnormally terminates with an
error message 'Failed to identify the eigenvalues' (which is very
unlikely to happen). However, in practice, appearance of this error
message indicates a divergence of iterations, which can be overcome by
fixing some of the Slater parameters.
The other two new parameters appearing on the first line of OUTGINE
are IORDER and ISAMPLE (columns 26-30 and 31-35, respectively).
If IORDER is non-zero, the eigenvector reordering with the
sampling algorithm will start only after n=IORDER iterations, and on
previous iterations, the eigenvectors will be sorted in the order of
increasing energy in each group of J values.
If ISAMPLE is non-zero, the list 6:26 PM 3/3/2013of eigenvectors after previous
iteration will be used as a sample for reordering eigenvectors after
iteration. Otherwise, reordering always uses a list of eigenvectors
from the zero-th (or, more precisely, IORDER-th) iteration as a sample
(coming from RCEINP or OUTGINE file, depending on the value of the
INVECT flag).

The same sampling algorithm is used in finding the level labels
(ascribing unique term names to levels), if the parameter CRIT1 is set
to a positive value (default is negative). This parameter is specified
in columns 61-65 of the line prior to the final '   -1' line in OUTGINE
for each parity (appeared since ver. 20 of RCE).

NOTE 1:

Changing the contents of columns 16-30 in the first line of OUTGINE for
the SECOND parity will have no effect. Only the very FIRST line of
OUTGINE is tested for this option.

NOTE 2:

For RCE iterations to work properly, the configuration labels used in
IN36 should be no longer than 5 characters and contain no spaces or
dashes. You can use primes, quotes, or asterisks to denote
configurations with differenct cores.

Description of RCEINP and RCEOUT file formats:
1) The line containing flags dealing with printing modes is the same as
in OUTGINE. Additional features produced by the same print flags:
a) Parameter IPRNA (if zero) makes all non-zero angular coefficients
for radial integrals (for each "experimental" level) to appear in the
output (OUTE), in much more readable form than it was before.
(Originally, this parameter caused printing of all the matrixes to be
diagonalized.)
b) Parameter IPRNV (if zero) makes all non-zero derivatives of
energy levels over the Slater parameters to be printed (originally, ALL
the derivatives matrixes were printed).
c) Parameter IPRNSQ=2 will suppress all the eigenvectors printout;
IPRNSQ=1 will suppress (as it were originally) only the output of
"eigenvector components squared".
2) The iteration control line (same as in OUTGINE).
3) The group of lines describing the eigenvectors resembles the format
of LEVELS1 file in LS-coupling: the energy; the asterisk (*) for
"unknown" level or empty space for "experimental" level; the calculated
(fitted) energy); J value; one or more groups describing the
eigenvector components: the amplitude with sign (multiplied by 100 and
rounded); configuration name and term name). Up to five components may
be present.
4) The group of lines with description of Slater parameters resembles
the format of OUTE file: parameter name, then the parameter flag,
than the value of parameter, then the maximum value for this parameter
(usually 0.000 - no max. value), then denominator (if DENOM flag is
set). For configuration-interaction parameters, names of configurations
are printed for reference. Do not change columns in which the numbers
are given (in the RCEINP file created during the first run of RCE).
The same structure (1-4) follows for the second parity (if present).

HARDWARE REQUIREMENTS (for PC users)
This package of COWAN programs needs IBM-compatible 386+387 or 486
computer with at least 1 Gb of contiguous free RAM. It cannot run on
386 machines without 387 co-processor. The programs may need up to
500 Mb of free hard disk space for normal execution, depending on
complexity and number of configurations computed simultaneously.

SOFTWARE REQUIREMENTS (for PC users)
This program package needs an operating system Windows XP, Vista,
Windows 7, or later. It was not tested on Windows 8. RCG and RCE
programs need at least 1 Gb of contiguous free RAM.


SUMMARY OF THE MAJOR FIXED BUGS
- Incorrect identification of Rydberg-series configurations in RCN2
led to appearance of non-zero Rd0 configuration-integrals in some
cases, e.g., when nsp3 and sp3n's sonfigurations are included in the
set, or when some core-excited configurations are present in addition
to a normal Rydberg series. See A. Kramida, Comput. Phys. Commun. 215,
47 (2017); Erratum to be published (2018).
- Multiple bug fixes made in the McGuinness version of RCN of 2004
have been incorporated in the present release of RCN. This solved poor
convergence problems in some cases, e.g., term-dependent calculations
in Ar I.
- Abnormal termination of RCG on complex calculations was caused in
earlier versions of the package by misalignment of some COMMON
block sizes, and a bug in the Salford Fortran linker resulted in
incorrect memory allocation for these COMMON blocks. This error has
been fixed by aligning all COMMON blocks having the same label to the
same size.
- Several dimension-checking blocks have been added to various
pieces of the RCG code.
- Added checking for the sufficiently large range of term sorting
keys in several places in RCG. If the sorting keys are insufficient,
a corresponding diagnostic ('Insufficient term sorting key') may
appear in the output. In such cases, RCG has to be re-compiled with
increased range of sorting keys (consult A. Kramida).
- Error that led to wrong JJ notations in the output files was
fixed in RCG and RCE.
- Error that led to creation of unexpected blank lines in the OUTGINE
file was fixed in RCG.
- Errors that sometimes led to appearance of control symbols in ING11
and OUTGINE files were fixed in RCN2 and RCG.
- RCN2 did not recognize the use of parameters with illegal range
(IABG=2, see the description of the input file for RCG), and if it was
necessary, the user needed to modify the ING11 manually. Now this bug is
fixed.
- In some cases, the original Cowan's RCG produced giant temporary
files (up to several Gbytes in size). Now they are reduced in size by a
factor of 10 to 100, speeding up the execution of RCG and RCE by a
factor of 2 to 5.
- In all programs, real*4 constants and functions (ABS, EXP, LOG
etc.) have been replaced with real*8 constants and corresponding
functions (DABS, DEXP, DLOG etc.). It is significant for some Fortran
compiles.
- Errors in RCG/RCE that led to wrong notation of levels in alternate
coupling schemes in RCE output files have been fixed.

RECOMPILING UNDER WINDOWS OS
The present executables have been compiled with the gfortran compiler 
bundled in the MinGW distribution package. I used v. 7.3.0 for the 
x86_64-w64 platform with Posix threads. More specifically, the package 
identifies itself as "gcc version 7.3.0 (x86_64-posix-seh-rev0, Built by 
MinGW-W64 project)." It can be downloaded from links in
https://gcc.gnu.org/wiki/GFortranBinaries#Windows. The Fortran source 
codes can be recompiled by using the batch file gfortran_O3.bat 
described in the ADDENDUM below.

COMPILING FOR OTHER OPERATING SYSTEMS
If you want to compile the Fortran texts for operation systems
different from Windows, you should eliminate all MSDOS-specific
modifications made in the RCG code for correct retrieval of the data
files SENIOR, FOR072, FOR073 and FOR074 on the disk. These modifications
are easy to find in first lines of RCG, where the files are opened. You
may simply change the 'MSDOS' flag to zero and it will do.
Another change that has to be made when compiling in any other
operating system is in the system-dependent function SECONDS() returning
the system clock timing in seconds. There are several variants suggested
for different systems. Choose the appropriate one and uncomment it,
and comment the PC-specific code (by placing 'C' in the first column
of all corresponding lines). This should be done for all programs of the
package. If you do not know the timing procedure for your specific
system, just make this function to always return zero.
The compilation itself is not sufficient for running RCG. You must
make the first run of RCG using a special input file (provided in the
package) to generate the CFP (coefficients of fractional parentage)
'decks' (files FOR072, FOR073 and FOR074). To do that, you should copy
the file ING11.CFP into the file ING11, delete the first 8 lines with
comments from it, and run RCG. The comment lines in ING11.CFP specify
the minimum array dimensions required for RCG to handle it.
Other programs of the package do not need any modifications or
special actions, except for the changes in function SECONDS() mentioned
above.

ADDENDUM: The list of files in the distribution package.

1) Root directory        - Documentation files (ASCII) copied from the
                           R. Cowan's directory at LANL on December 1994,
                           except when otherwise indicated.
Contents:
        RCN_DOC.txt        - manual for RCN/RCN2 and HF8 programs;
        RCG_DOC.txt        - manual for RCG program;
        RCE_DOC.txt        - manual for RCE program;
        README.Cowan.htm   - R. Cowan's write-up for his version of the
                             package of December 1999;
        README_Kramida.txt - this file;
        license.txt        - NIST license and warranty information.

2) Code directory - executable files and input files common for all
                     calculations.
Contents:

a) Main executables:
        RCN36K.EXE - Program to compute single-config. wavefunctions
        RCN2K.EXE  - Program to compute single-config. and CI radial
                    integral and transition matrix elements
        RCG11K.EXE - Main program to compute spectra, cross-sections etc.
        RCE20K.EXE - The least squares fitting program.

b) Main batch files for running the main executables:
        RCN.BAT       - batch file for running the RCN program;
        RCN2.BAT      - batch file for running the RCN2 program;
        RCG.BAT       - batch file for running the RCG program;
        RCE.BAT       - batch file for running the RCE program;
        make_cfp.bat  - batch file for creating the binary CFP decks.

c) Input files common for all calculations:
        FOR072        | binary files with CFP decks used by RCG
        FOR073        | (not included in the distribution package;
        FOR074        | these files are created during installation);

        SENIOR        - ASCII file used by RCG;

d) gfortran executable files and binary libraries:
        gfortran_O3.bat    - Batch file for compiling the source codes with
                       the GNU gfortran compiler. Usage:
                       gfortran_O3 <source_file_name>
                       (The file name must be given without extension,
                       which must be .F)
        libquadmath-0.dll  - gfortran binary library.
        libgcc_s_seh-1.dll - gfortran binary library.

e) Auxiliary files (rarely needed):
        ING11.CFP     - the complete CFP input deck for RCG.

f) Utility files:

        f.1) Perl utilities:
                       Each of the .bat files is accompanied by a
                       corresponding .pl Perl source file. You can
                       run the .pl files separately by

                       perl <file_name.pl> [parameters]

                       or in debug mode

                       perl -d <file_name.pl> [parameters]

                       Two Perl programs, print_out.pl and conv_out.pl, 
                       require an additional option to be specified:

                       -I<path_to_COWAN\CODE>

                       Otherwise, Perl will not find the included
                       file located in the COWAN\CODE directory.

                       You can modify any of the .pl files as needed.
                       However, if you find any bugs or introduce
                       new useful options, please notify me, so that
                       your changes and bug fixes could be
                       incorporated in future releases.

        ADD_EAV.BAT   - Perl program used to add a given energy shift
                       to all Eav parameters in the ING11 file.
                       Usage:

                       add_eav <energy_shift_in_kK>

                       The same action can be achieved by including
                       the rescale 'card' in ING11, i.e. by adding
                       at the top of ING11 a line with zero in
                       position 5 and the value of the energy shift
                       in positions 21 through 30.

        update11.bat  - A batch file to be executed after the LSF fitting
                       with RCE to transfer the fitted parameters into ING11,
                       run RCG with those fitted parameters, and run
                       printout.bat (see below). If the calculation includes 
                       both parities, this command must include "2" as a 
                       parameter:

                       update11 2

                       This is the easiest and safest way to produce 
                       Excel-readable list of energy levels and LSF 
                       parameters. It must also be executed before running 
                       the conv_out utility code (see below).

        S11.BAT       - Perl program for substituting the parameters
                       from the RCEOUT file (produced by RCE
                       least-squares fitting program) back into the
                       ING11 file (input for RCG program);
                       Usage:

                       s11

                       (with no parameters - for calculations with
                       one parity set or for substituting the
                       fitted parameters into ING11 for the first
                       parity set only);

                       s11 2

                       (for substituting the fitted parameters into
                       ING11 for two parity sets).

        printout.bat  - Perl program to print out the level list in
                       a form similar to LEVELS2, but with
                       configuration and term labels including
                       proper parent terms in the order of summation
                       of shells. This works more-or-less well with
                       LS and JJ coupling only. The eigenvector
                       compositions with parentage are read from the
                       OUTG11 file, while the experimental energies
                       are substituted from the RCEOUT file. Thus,
                       to produce correct results, you must first
                       substitute the fitted parameters from RCEOUT
                       into ING11 (e.g., by running s11 2) and
                       re-run RCG. The program also produces a file
                       with fitted and HF parameter values and ratios
                       LSF/HF. Usage:

                       printout <levels_out_name> [params_out_name]

                       where <levels_out_name> is the name of the
                       output levels file (typically, levels.txt),
                       and [params_out_file] is the optional name
                       of the parameters output file (if omitted,
                       the default is params.txt). The output files
                       are ASCII tab-delimited files best viewed
                       with Excel or any other spreadsheet program.

                       This program takes into account the possible
                       reordering of subshells (see description of
                       reorder_ing11.bat below). For it to work,
                       the files in36, OUTG11, OUTE, RCEINP.HF, and
                       LEVELS1 must be present in the current
                       directory. The scaling factors for the HF
                       parameters are read from OUTG11. Therefore,
                       for the ratios LSF/HF to be correct, the
                       first run of RCE, which produces RCEINP.HF,
                       must be made with the sequence of runs
                       RCN/RCN2/RCG/RCE. The files OUTE and RCEINP.HF
                       are skipped if you specify an optional /noRCE
                       parameter in the command line.

        conv_out.bat  - Perl program to convert both level and line
                       lists from the output of RCG and RCE. Input
                       parameters:
                       1) Output levels file name (e.g., 'levels.txt').
                          The format of this particular file is 
                          slightly different from a similar file created 
                          by the printout code. Namely, it contains 
                          several additional columns. The most important 
                          one is the unique level number, which is used 
                          to label the lower and upper levels in the 
                          simultaneously created transitions file.
                       2) Output transitions file name (e.g., tr.txt)
                       3) RCE or noRCE
                          (for including or omitting RCE data such
                          as experimental energy levels)
                       4) LTE or DR or AUGER or CM or
                          BF:<time_of_flight_ns>

                          This specifies how to scale the line
                          intensities. The most popular use is LTE
                          (local thermodynamic equilibrium, a.k.a.
                          Boltzmann). DR stands for 'Dielectronic
                          Recombination'. This option, as well as
                          AUGER requires autoionization rates to be
                          available in two separate files, one for
                          the first parity set and one for the second.
                          If this option is specified, the names of
                          these files must be given at the end of the
                          parameter list (parameters 7 and 8).

                          CM is a special 'Cascade Matrix' option
                          in calculations of Auger electron spectra.

                          BF stands for 'Beam-foil'. For this option,
                          an additional parameter, time of flight
                          in nanoseconds, must be given together
                          with 'BF:'. It means the time required for
                          the emitting ions to cross the field of
                          view of the spectrometer. The intensity
                          is estimated by integrating the line
                          intensity from each of the decaying levels
                          within this 'time of flight' period.
                       5) Effective excitation temperature in 1000
                          cm-1. This parameter is required with all
                          options, but is used only in LTE.
                       6) Max. scaled intensity
                          A number meaning the ratio of max and min
                          intensities of the lines to be included in
                          the transition list. The min. intensity is
                          assigned the value of 1, and the rest are
                          scaled accordingly. Lines with intensity
                          smaller than 0.5 are omitted. All
                          intensities are rounded to an integer
                          number.
                       7,8) Optional parameters: Autoionization
                          rate files for the first and second parity
                          configuration sets. These files are
                          generated by a special Perl utility code
                          read_aa.bat. It is rater difficult to
                          give a proper description of this code.
                          If you are able to understand the Perl
                          source code (which includes some comments
                          and examples of usage), you can use it.
                          Otherwise, email alexander.kramida @
                          nist.gov for advice or do not use the DR
                          and AUGER options.
                        9) Optional parameter: /IDEN
                          It is used to generate input files 
                          ENLEV.DAT and TRANS.DAT for the visual
                          line identification code IDEN2 
                          (V. I. Azarov, A. Kramida, and M. Ya. 
                          Vokhmentsev, Comput. Phys. Commun. 225, 
                          149–153 (2018); 
                          DOI:10.1016/j.cpc.2017.12.012). This 
                          will work only if the observed line list
                          DLV.DAT and an auxiliary file numset.dat
                          are present in the working directory
                          (for explanation of those files, see the
                          above article on the IDEN2 code).

                       In the transitions output file, with the LTE,
                       BF, and DR options, the wavelengths are given
                       in angstroms. Those between 2000 and 20000
                       angstroms (or, more precisely, those
                       corresponding to wavenumbers between 5000 and
                       50000 cm-1) are given in standard air,
                       otherwise in vacuum. With the AUGER options,
                       the transition energy in eV is printed
                       instead.

                       For each line, the configuration and
                       LS term label are printed for both levels,
                       as well as the Lande factors, radiative
                       lifetimes, gA and A values, the value of the
                       cancellation factor, and, if autoionization
                       data are included, the corresponding
                       autoionization rates and branching fractions.

                       Admittedly, this Perl version of the old
                       Borland Pascal code conv_out.exe is much less
                       developed and lacks many useful options.
                       However, it has some advantages, such as no
                       limitation on the number of transitions and
                       energy levels, and portability. In the future,
                       it may completely replace the old Pascal
                       version.

        read_aa.bat
                     - Perl program to process the output of RCG
                       for autoionization calculations (see above).

        reorder_ing11.bat
                     - Perl program for reordering the subshells in
                       the ING11 file. This is needed if you want to
                       change the default summation order of the
                       subshells. For example, the configuration
                       1s.2s.2p of Li-like spectra is best described
                       by first combining the 2s and 2p subshells to
                       form intermediate LS terms 3P* or 1P*, and
                       then combining this intermediate term with
                       the 2S term of the 1s subshell. The default
                       order of shell summation in RCG is left to
                       right, i.e., 1s is combined with 2s to
                       produce 3S of 1S, then this is combined with
                       the remaining 2p shell to produce the final
                       4P* or 2P* term. Usage:

                       reorder_ing11 <output_filename> <new_order>

                       where <output_filename> is the name of the
                       new ING11 file (the original ING11 must be
                       present in the current directory; it is used
                       as input), and <new_order> is a sequence of
                       digits, e.g., 231 to specify that shells 2
                       and 3 will go first in the new ING11,
                       followed by shell 3. You can recognize the
                       shells in the top section of ING11 and count
                       their numbers from left to right, starting
                       with 1. The program rearranges not only the
                       shell designations, but also the corresponding
                       Slater and CI parameters.

        scale_param.bat
                     - Perl program for scaling the parameter values
                       in the RCEINP file by given factor(s).
                       Parameters:
                       1) input_file_name
                       Name of the input file (normally, a saved copy
                       of the RCEINP file).
                       2) output_file_name
                       Name of the new RCEINP file to produce.
                       3) One or more of the following tokens
                       separated by space:
                       F=<F_factor>
                       G=<G_factor>
                       CI=<CI_factor>
                       GROUP<number>=<factor>
                       where the factors given are used to multiply
                       all parameters from the corresponding group
                       of the input file. The GROUP number is the
                       absolute value of the numerical flag used to
                       group parameters in RCEINP, e.g., 51, 52, etc.
                       4) conf_set
                       (1 for the first parity set, 2 for the second
                       parity set, or 12 for both).
                       5) FIX
                       An optional parameter. If specified, all Slater
                       parameters that were 'free' in the input file
                       will be 'fixed' in the output file for the
                       given parity set(s).

        subst_exp_levs.bat
                     - Perl program to substitute

                       Usage:

                       subst_exp_levs <input file> <output_file> [params=flags|noflags|Eav|F|G|Z|CI|none]

                       This program is used to substitute experimental
                       level values from one RCEINP (or RCEOUT) file
                       into another, possibly created with a different
                       set of configurations. Typically, you would
                       need it if you made an LSF fitting with a
                       limited set of configurations and then
                       discovered you need more configurations to
                       include. Then, after you run RCN/RCN2/RCG/RCE
                       with an extended configuration set, instead of
                       manually editing again your new RCEINP file,
                       you can use this program to transfer the level
                       values from your previous calculations to the
                       new RCEINP file, which is the <output_file> in
                       the parameter list.

                       Optionally, you can specify a certain action
                       for the parameters from your old calculation.

                       params=flags will set only the flags for the
                       same parameters in the new RCEINP (i.e., make
                       them fixed or free or grouped in the same way
                       as they were in your old RCEINP).

                       params=noflags will transfer only the values
                       of the parameters but not the flags.

                       params=Eav or params=F or params=G or params=Z
                       or params=CI will transfer only the values and
                       flags for the corresponding parameter group
                       (Z stands for ZETA). Both values and flags
                       will be transferred.

                       params=none will skip all parameter values and
                       flags (leave them unchanged in the new RCEINP
                       file).

        subst_params_ing11.bat
                     - Perl program to substitute parameter values
                       from one ING11 file into another one. The
                       substitution is done line by line. If the
                       configuration name(s) on the line read from the
                       input file matches those in the output file, the
                       parameter line(s) for these configurations are
                       substituted in place of similar lines in the
                       output file. The input ING11 file is assumed to
                       have both parities, while the output ING11 is
                       supposed to have only one parity set. The main
                       purpose is to substitute fitted parameters to
                       set up ING11 for autoionization calculations.

                       Usage:

                       subst_params_ing11 <input_ING11> <output_ING11> <conf_set>

                       where <conf_set> can be 1 or 2, for conf. set #1
                       or #2 in the <input_ING11>.

        conf_list.bat - Perl program to print out the lists of
                       configurations included in the first and second
                       parity sets of IN36.

        vacair.pl     - An include file for conv_out.pl (and .bat)
                       containing functions for conversion of
                       wavelengths from vacuum to standard air and
                       vice versa. The function names are Lvac() and
                       Lair(), and their arguments are the air or
                       vacuum wavelengths in angstroms, respectively.


        f.2) Old Borland Pascal utilities:
        These executables no longer work directly under Windows 7 or
        any other 64-bit version of Windows. However, they still can be
        used in conjunction with the freeware DosBox utility.

        HCI.BAT       - batch file for running HIGH_CI.EXE program
                       (to determine the list of highly interacting
                        configurations from data in file ING11);
        HIGH_CI.EXE   - program used by HCI.BAT file;
        CO.BAT        - batch file for conversion of the calculated
                       transition probabilities and energy levels
                       from files OUTG11 and LEVELS2 (output files
                       of the RCG and RCE programs) into the input
                       files for V. Azarov's IDEN program (for lines
                       identification);
        CONV_OUT_old.EXE
                     - exe-file used by CO.BAT conversion program;
        CONV_OUT.HLP  - help on usage of CO.BAT conversion program;
        CONV_OUT.PAR  - parameter file used by CO.BAT conversion program;
        FINDDIFF.EXE  - program used for comparison of level lists
                       for preparation of input file for MTRANS1.EXE)
        MTRANS1.EXE   - program for replacement of calculated energy
                       levels and transition probabilities in the input
                       files for V. Azarov's IDEN program with those from
                       subsequent runs of RCE/RCG (keeping all the
                       identified lines and levels intact);

3) FOR directory    - Fortran 77 source files of the Cowan programs.
        RCN36K.F
        RCN2K.F
        RCG11K.F
        RCGPAR.F      - Include file for RCG11K.F
        RCE20K.F
        RCEBPAR.FOR   - Include file for RCE20K.F.


4) WORK directory   - Sample directory for Cowan calculation projects.
                       Includes sample input files for simple Sn VIII
                       calculations with a single configuration in each
                       parity set.

        IN36          - sample input file for RCN (for two configurations
                       of Sn VIII);
        IN2           - sample input file for RCN2;


                                  April 21, 2021

                                  Alexander Kramida,
                                  tel. +1 (301) 975-8074
                                  e-mail alexander.kramida @ nist.gov

LICENSE AND WARRANTY INFORMATION
