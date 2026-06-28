{$A+,B-,E-,F+,G-,I-,N+,O-,R-,S-,V-,X-}
{$M 16384,0,655360}
program Convert_RCG_output_to_IDEN;
{     Program Conv_Out is used for conversion of the RCG and RCE
output files into the IDEN input files TRANS.DAT and ENLEV.DAT.

Besides it, this program creates reduced linelist, sorted by line wavelengths
(this option did not work properly in original RCG).

All input parameters needed for the CONV_OUT program must be specified
in the parameter file named CONV_OUT.PAR. It must be in the current directory
when you run CONV_OUT. Copy the sample CONV_OUT.PAR to your working directory
and modify it as you need. The sample file contains description of its format
inside. All comments written in the parameter file must start from column 40
or more.

What the CONV_OUT program does:

    - Reads the parameter file CONV_OUT.PAR;

    - Reads the OUTG11 file created by RCG mod 11;
    - Reads the RCEOUT file created by RCE mod 20k;

    CAUTION: RCE and RCG MUST be ran with the same
    radial parameters (RCE - with zero number of iterations,
    just after the RCG, or RCG with ALL parameters taken from the
    last iteration of RCE), to get the most close eigenvalues and
    eigenvectors. RCG should be ran so that OUTG11 contains transition
    list sorted by energy levels (of first and/or second parity);

    - Analyzes OUTG11 (to get the maximum line intensity and minimum
      radiative lifetime) and RCEOUT. The program reads the eigenvector
      amplitudes from OUTG11 and from RCEOUT. The level energies and
      compositions in both files are compared with each other to get the
      closest match.
         Normally, if you have any warning or abnormal termination
      during the Analize stage, you should make sure that RCE and RCG
      have been ran with the same radial parameters (see above).

    - Reduces the line list, leaving alive only transitions
      that are more intense than given threshold value;
      this line list will be wavelength-sorted. Transition probabilities
      are printed as gA*Scale/gAmax, where the Scale factor is given
      in the parameter file.
      You may use two ways of specifying the minimum line intensity,
      starting from which, the line will appear in the linelist. You may
      give the minimum absolute value of transition probability in sec-1,
      and/or the minimum scaled relative intensity value.

    - in this "reduced linelist", you can artificially increase
      intensities of lines arising down from metastable levels by
      some factor - to make some account for the
      rise of their population; you should idnicate explicitly,
      what levels should be treated as metastable. For that,
      specify (in the parameter file) the minimum level lifetime to treat
      the level as metastable. If any level has longer lifetime than this
      minimum value, it will be treated as metastable: the intensities of
      all lines arising down from it will be multiplied by a factor
      indicated in the parameter file. If this factor is <= 1, this
      option does not work. These scaled lines are written also to
      a separate file named "META.", as well as the list of metastable
      levels.

    - Converts the "reduced linelist" to the files ENLEV.DAT and
      TRANS.DAT for the V.Azarov's IDEN program.
      You can remove the levels which are not connected with others by
      any transition, reducing the total number of levels.

    - Creates (by your will) the energy levels file sorted in a special
      way named CONF-J-E. All levels are first divided into groups by
      configuration names. Each group is divided into sub-groups with
      the same J. Then each of these sub-groups is sorted by level
      energies.
      This file is needed if you are upgrading your COWAN calculations
      using newly attained experimental levels data. After "improved"
      calculations, many levels will change designations, and their
      energy order will differ from what you started from. But your
      old ENLEV.DAT and TRANS.DAT still contains many identifications
      that you have made to get the new levels data. Not to loose this
      work, you may transfer these identifications to the new ENLEV.DAT
      and TRANS.DAT. But for this, you must install the correspondence
      between the level numbers in old ENLEV.DAT and new ENLEV.DAT.
      The CONF-J-E sorted file will help you to do it. You should also
      have the CONF-J-E sorted file produced from the old ENLEV.DAT file,
      and bring the two files together, new data in the left column and
      old data in the right column. The separate program FINDDIFF checks
      this combined file, making the list of lines having differences in
      level notations. This list usually is not very long, and you can
      work through it, replacing halves of lines so that in every line the
      left half and the right half contain the same electron state.
      After that, the program MTRANS1 can produce the combined file
      TRANS.DAT containing all old identifications but new calculated
      levels and transition probabilities.


                                           A.E.Kramida,
                                           July 3, 1993                    }
uses Objects,VacAir;
const
  MaxJval  =15;  { Max. number of different J values in each parity }
  MaxAmplComps=5; { Max. number of amplitude components for level in RCEOUT}
  MaxCollectionSize=16380;

type
  real=single;
  str15 = string[15];

  PAmpl=^TAmpl;
  TAmpl=object(TObject)
    Ampl: Integer;
    ConfNo: Integer;
    TermRef: PAmpl;
    TermName: string[5];
  end;

  StrMaxJval=string[MaxJval];
  PLevel = ^TLevel;
  TLevel = object(TObject)
    N: integer;
    E,ET: real;
    Lft: extended;
    AA: PCollection;   { gA values for transitions to other levels }
    Ampls: PCollection; { Map of amplitudes for composition }
    gList: ^StrMaxJval;
    g: byte;
    Fixed: boolean;
    constructor Init;
    destructor Done; virtual;
    function GetBestE: real;
    {function Exp: boolean;}
  end;

  PConf =^TConf;
  TConf = object(TObject)
    Name: string[5];
  end;

  PLine=^TLine;

  PTrans=^TTrans;
  TTrans=object(TObject)
    gA: real;        { Will be converted to intensity scale }
    gA0: real;       { Will stay equal to original Cowan's gA }
    gf: real;       { Will stay equal to original Cowan's gf }
    Q: boolean;
		cf: real;       { Cancellation factor }
    FromLevel,
    ToLevel: PLevel;
    ObsLine: PLine;
  end;

  TLine=object(TObject)
    Int: integer;
    Wn: double;
{    Comments: string[40];}
  end;

  PNumSetLine=^TNumSetLine;
  TNumSetLine=object(TObject)
    Start: integer;
    dWL: single;
  end;

  PRange=^TRange;
  TRange=object(TObject)
    FromWL, ToWL: double;
    gAScale: double;
    gAmax: double;
    Absorption: boolean;
    Teff: double;
  end;

procedure CollectionSizeError(S: string);
begin
  writeln('Max. collection size exceeded for ',S);
  halt(1);
end;

procedure Trim(var S: string);
begin
  while (S[1]=' ') and (S<>'') do delete(S,1,1);
  while (S<>'') and (S[length(S)]=' ') do dec(S[0])
end;

constructor TLevel.Init;
begin
  if not TObject.Init then
    Fail;
  AA:=New(PCollection,Init(1,1));
  if AA=nil then
    Fail;
  Ampls:=New(PCollection,Init(1,1));
  if Ampls=nil then
    Fail;
  gList:=nil;
  Fixed:=false
end;

destructor TLevel.Done;
begin
  Dispose(AA,Done);
  Dispose(Ampls,Done);
  if gList<>nil then
    FreeMem(gList,SizeOf(gList^))
end;

function TLevel.GetBestE: real;
begin
  {if (E<0) and (ET>0) or (E>0) and (ET<0) then}
  if not Fixed then
    GetBestE:=ET
  else
    GetBestE:=E
end;

(*
function TLevel.Exp: boolean;
begin
  Exp:={(E>=0) and (ET>=0) or (E<=0) and (ET<=0)}fixed
end;
*)

var
  fRCGout,f2,f3,f4,f5,fRCEOUT,fm,fNUMSET,fDLV: text;
  S,Sx: string;
  S1,S2,S3,CName: string[15];
  Res,Nc,i1,Ipar,gprev,Ir,NJv,IL,IL1,IL2,NRanges: integer;
  NJval: array[1..2] of integer;     { Number of J-values in each parity }
  NLevJ: array[1..MaxJVal,1..2] of integer; { Number of levels in each J-group}
  gA,AgrMax: real;
  MS,Start,Stop,Stop1,UseNames,MakeSLFile,
  RemoveUnconnected,TabDelimitedOutput: boolean;
  k,n,n1,p1,p2,i,j,j1,jm,Irev,Ij,Nlev,Nu,NL,NLms,NTms,Nms,ND,Nfound,
  BestMatch: word;
  lambda,Imin,Lmin,Lms,Fms,dEexp,dEth,Eu,El,E1,E2,Br,Aj,
  MinTermSplit,Athresh,dE,MinAirWn,MaxAirWn: real;
  A1,A2: integer;
  L: extended;
  Dist,MinDist: LongInt;
  Pu,Pl,PL1,PGround: PLevel;
  PR: PRange;
  ch: char;
  CT: array[1..1000] of PConf;
  PC: PConf;
  PAm,PAm1: PAmpl;
  Tr: PTrans;
  Second,Header,Paint,SamePar: boolean;
  Kprint: byte;
  RCELevels,RCGLevels,Trans,NumSet,Lines,ResTrans: TCollection;
  LSBasis: array[1..2] of TCollection;
  Ranges: TCollection;
  Emax: real;
const Absorp: array[boolean] of string[10]= ('Emission','Absorption');
label 1;

procedure Error(S: string);
begin
  writeln;
  writeln(S);
  halt(1)
end;

procedure ReadParams;
var
  F: text;
  S: string;
procedure ReadS;
var p: byte;
begin
  readln(f,S);
  if IOResult<>0 then
    Error('Error reading parameters file.'^G);
  p:=pos(';',S);
  if p=0 then
    p:=40;
  S:=copy(S,1,p-1);
  Trim(S)
end;
begin   { ReadParams }
  assign(f,'conv_out.par');
  reset(f);
  if IOResult<>0 then begin
    writeln('Can not open parameters file CONV_OUT.PAR.'^G);
    halt(1)
  end;

  readS;
  if S='' then
    Error('Transitions file name missing in parameters file.'^G);
  assign(fRCGout,S);
  reset(fRCGout);
  if IOResult<>0 then
    Error('Can not open transitions file.'^G);

  readS;
  if S='' then
    Error('RCEOUT file name missing in parameters file.'^G);
  assign(fRCEOUT,S);
  reset(fRCEOUT);
  if IOResult<>0 then
    Error('Can not open RCE vector output file.'^G);

  readS;
  if S='' then
    Error('Reduced linelist file name missing in parameters file.'^G);
  assign(f2,S);
  rewrite(f2);
  if IOResult<>0 then
    Error('Can not create reduced linelist file '+S+^G);

  ReadS;
  val(S,Imin,Res);
  if Res<>0 then
    Error('Syntax error in parameters file (min. intensity line).'^G);

  ReadS;
  val(S,Emax,Res);
  if Res<>0 then
    Error('Syntax error in parameters file (max. level-energy line).'^G);
  Emax:=Emax*1e-3;

  ReadS;
  val(S,NRanges,Res);
  if Res<>0 then
    Error('Syntax error in parameters file (No. of ranges line).'^G);
  if NRanges<=0 then
    Error('Illegal value for number of ranges in parameter file.'^G);

  for i:=1 to NRanges do begin
    if MaxAvail<SizeOf(TRange)+1000 then
      Error('Insufficient memory.'^G);
    PR:=New(PRange,Init);
    if PR=nil then
      Error('Low memory.'^G);
    ReadS;
    val(S,PR^.FromWL,Res);
    if Res<>0 then
      Error('Syntax error in parameters file (min. wavelength line).'^G);
    ReadS;
    val(S,PR^.ToWL,Res);
    if Res<>0 then
      Error('Syntax error in parameters file (max. wavelength line).'^G);
    ReadS;
    val(S,PR^.gAScale,Res);
    if Res<>0 then
      Error('Syntax error in parameters file (intensity scale line).'^G);
    PR^.gAmax:=0;
    ReadS;
    if not (upcase(S[1]) in ['Y','N']) then
      Error('Syntax error in parameters file (absorption-emission line).'^G);
    PR^.Absorption:=upcase(S[1])='Y';
    readS;
    val(S,PR^.Teff,Res);
    if Res<>0 then
      Error('Syntax error in parameters file (eff. temperature line).'^G);
    PR^.Teff:=PR^.Teff*1e-3;

    Ranges.Insert(PR)
  end;

  readS;
  ch:=upcase(S[1]);
  if ch<>'Y' then begin
    Fms:=1;
    Lms:=9e31;
    ReadS;
    ReadS
  end else begin
    ReadS;
    val(S,Lms,Res);
    if Res<>0 then
      Error('Syntax error in parameters file (min lifetime line).'^G);
    ReadS;
    val(S,Fms,Res);
    if Res<>0 then
      Error('Syntax error in parameters file (metastables scaling line).'^G)
  end;

  readS;
  val(S,MinTermSplit,Res);
  if Res<>0 then
    Error('Syntax error in parameters file (min. term splitting line).'^G);
  MinTermSplit:=MinTermSplit*1e-3;

{  ReadS;
  val(S,KPrint,Res);
  if Res<>0 then
    Error('Syntax error in parameters file (print sorting line).'^G);}

  ReadS;
  val(S,i,Res);
  if Res<>0 then
    Error('Syntax error in parameters file (configurations naming line).'^G);
  UseNames:=i<>0;

  ReadS;
  val(S,dEexp,Res);
  if Res<>0 then
    Error('Syntax error in parameters file (exp. levels uncertainty line).'^G);

  ReadS;
  val(S,dEth,Res);
  if Res<>0 then
    Error('Syntax error in parameters file (th. levels uncertainty line).'^G);

  readS;
  ch:=upcase(S[1]);
  if not (ch in ['Y','N']) then
    Error('Syntax error in parameters file (remove unconnected levels line).'^G);
  RemoveUnconnected:=ch='Y';

  readS;
  MakeSLFile:=S<>'';
  if MakeSLFile then begin
    assign(f5,S);
    rewrite(f5);
    if IOResult<>0 then
      Error('Can not create CONF-J-E sorted levels file '+S+^G)
  end;

  readS;
  ch:=upcase(S[1]);
  if not (ch in ['Y','N']) then
    Error('Syntax error in parameters file (paint matching lines line).'^G);
  Paint:=ch='Y';

  ReadS;
  if Paint then begin
    val(S,Athresh,Res);
    if Res<>0 then
      Error('Syntax error in parameters file (assignment threshold line).'^G);
  end;

  ReadS;
  val(S,MinAirWn,Res);
  if Res<>0 then
    Error('Syntax error in parameters file (min. air wavenumber).'^G);

  ReadS;
  val(S,MaxAirWn,Res);
  if Res<>0 then
    Error('Syntax error in parameters file (max. air wavenumber).'^G);

  ReadS;
  ch:=upcase(S[1]);
  if not (ch in ['Y','N']) then
    Error('Syntax error in parameters file (Tab-Delimited Output line).'^G);
  TabDelimitedOutput := ch='Y';

  close(f)
end;  { ReadParams }

function GetConfName(S1: string): string;
var i,delta: byte;
begin
  while (S1<>'') and (S1[length(S1)]=' ') do dec(S1[0]);
  if length(S1)>5 then begin
    repeat
      i:=pos(' ',S1);
      if i<>0 then delete(S1,i,1)
    until (i=0) or (S1='');
    S1:=Copy(S1,1,5)
  end;
  if S1<>'' then begin
    delta:=ord('a')-ord('A');
    for i:=1 to length(S1) do
      if S1[i] in ['A'..'Z'] then
        inc(S1[i],delta);
    while length(S1)<5 do S1:=S1+' '
  end;
  GetConfName:=S1
end;

function GetConfNo(var S1: string): integer;
var i,j: integer;
begin
  j:=-1;
  for i:=1 to 1000 do
    if (CT[i]<>nil) and (CT[i]^.Name=S1) then begin
      j:=i;
      break;
    end;
  GetConfNo:=j
end;

function GetG(var Js: string; var Res: Integer): byte;
var R: real;
begin
  Trim(Js);
  if Js[1]='.' then insert('0',Js,1);
  val(Js,R,Res);
  GetG:=round(2*R+1)
end;

function GetG1(var Js: string; var Res: Integer): byte;
var R,p: integer;
begin
  Trim(Js);
  p:=pos('/',Js);
  if p<>0 then
    Js[0]:=chr(p-1);
  val(Js,R,Res);
  if p=0 then
    GetG1:=2*R+1
  else
    GetG1:=R+1
end;

function GetLevel(var S: string; var Pu: PLevel): boolean;
var i,j: word;
  p,j1: byte;
  ini: boolean;
  CC,Js: string[15];
  R: real;
  Nc: integer;
  PAm: PAmpl;
{procedure Est(P: PLevel); far;
begin
  if (J1=P^.g) and (j=P^.ConfNo) and (CC=P^.Term) and
    ((Pu=nil) or (abs(Eu-P^.ET)<abs(Eu-Pu^.ET))) then
        Pu:=P
end;}

procedure Est(P: PLevel); far;
var Pam: PAmpl; Nc: integer;
begin
  if (J1<>P^.g) then
    Exit;
  PAm:=P^.Ampls^.At(0);
  Nc:=PAm^.ConfNo;
  if PAm=nil then
    Exit;
  if not SamePar then
    if ini then begin
      if (not Second) and (Nc>500) then
        Exit;
      if Second and (Nc<=500) then
        Exit;
    end else begin
      if (not Second) and (Nc<=500) then
        Exit;
      if Second and (Nc>500) then
        Exit;
    end;
  if (Pu=nil) then begin
    Pu:=P;
    Exit
  end;
  if (abs(Eu-P^.ET)<=abs(Eu-Pu^.ET)) and (Nc=j)
  and (CC=PAm^.TermName) then
    Pu:=P
end;

begin { GetLevel }
  GetLevel:=false;
  p:=pos('* * *',S);
  ini:=p<>0;
  if ini then begin
    { Set correct value for the parity of initial level }
{          E      J     LEVEL        DELTA E  LAMBDA(A)  LOG GF GA(SEC-1) CF/BR

    * * *      88.452  2.0   1 (4S) 5S      NeV   2s2p3        * * *

   1      0.196  1.0  1 (3P) 3P       88.256   1133.067  -5.564 1.419E+04 -0.76
   2    413.443  1.0  2 (3P) 3P      324.991    307.700  -5.239 4.066E+05  0.44
   3    635.928  1.0  3 (2P) 1P      547.476    182.656  -8.377 8.386E+02  0.22}
    S1:=GetConfName(copy(S,p+46,7));
    j:=GetConfNo(S1);    { j is the conf.no in CT }
    Second:=j>500;    { set initial parity to Second if j>500 }

    Js:=copy(S,p+18,4);      { get J of the level }

    CC:=copy(S,p+28,2)+copy(S,p+32,3); { Term name }
    S1:=copy(S,p+5,12)                 { Energy }
  end else begin
    p:=pos('.',S);
    S1:=copy(S,p+10,2);
    Trim(S1);
    val(S1,j,Res); { get conf. number }

    CC:=copy(S,p+14,2)+copy(S,p+18,3);
    if (not Second) and (not SamePar) then
     { search for j-th conf. among second parity }
      j:=1001-j;
    Js:=copy(S,p+5,4);
    S1:=copy(S,p-6,10)
  end;
  Trim(CC);

  j1:=GetG(Js,Res);

  Trim(S1);
  if Res=0 then
    val(S1,Eu,Res);

  if Res<>0 then
    Error('Format error in transitions file.'^G);

  Pu:=nil;
  RCGLevels.ForEach(@Est);

  GetLevel:=Pu<>nil;

  if (Pu=nil) or (abs(Eu-Pu^.ET)>5e-1) then begin
    writeln;
    writeln('Error: level ',Eu:12:3,' not found in the list.');
    if Pu<>nil then
      writeln('Best close E is ', Pu^.ET:12:3);
    halt(1)
  end;

  PAm:=Pu^.Ampls^.At(0);
  Nc:=PAm^.ConfNo;
  if (Nc<>j) or (CC<>PAm^.TermName) then begin
    Writeln('Warning: level name mismatch in RCG: E=',Eu:12:3,
            '; J=',(Pu^.g-1)/2:4:1);
    writeln('Nconf = ',j, '; Ref. Nconf = ',Nc,'; Term = ', CC, '; Ref. term = ',
      PAm^.TermName);
  end;
  if abs(Eu-Pu^.ET)>0.01 then
    Writeln('Warning: level energy mismatch in RCG: E=',Eu:12:3,
            ' J=',(Pu^.g-1)/2:4:1);
end;

procedure ReadCTable;
var p1,j: byte; i,i1,i2,i3: integer;
procedure ReplaceCN;
procedure ReplaceCN1(PL: PLevel); far;
var PAm: PAmpl; k: integer;
begin
  for k:=0 to PL^.Ampls^.Count-1 do begin
    PAm:=PL^.Ampls^.At(k);
    if (PAm<>nil) and (PAm^.ConfNo=i2) then
      PAm^.ConfNo:=0
  end
end;
procedure ReplaceCN2(PL: PLevel); far;
var PAm: PAmpl; k: integer;
begin
  for k:=0 to PL^.Ampls^.Count-1 do begin
    PAm:=PL^.Ampls^.At(k);
    if PAm<>nil then
      if PAm^.ConfNo=i1 then
        PAm^.ConfNo:=i2
      else if PAm^.ConfNo=0 then
        PAm^.ConfNo:=i1
  end
end;
procedure ReplaceCN3(PAm: PAmpl); far;
begin
  if PAm<>nil then
    if PAm^.ConfNo=i1 then
      PAm^.ConfNo:=i2
    else if PAm^.ConfNo=i2 then
      PAm^.ConfNo:=i1
end;
var i: byte;
begin { ReplaceCN }
  if i1<>i2 then begin
    RCGLevels.ForEach(@ReplaceCN1);
    RCGLevels.ForEach(@ReplaceCN2);
    for i:=1 to 2 do
      LSBasis[i].ForEach(@ReplaceCN3);
    PC:=CT[i2];
    CT[i2]:=CT[i1];
    CT[i1]:=PC
  end
end;

begin  { ReadCTable }
{
         ELEC DIP SPECTRUM     (ENERGIES IN UNITS OF  1000.0 CM-1 =    0.12 EV)

   1      NeV   2p2            ---      NeV   2s2p3
   2      NeV   2p4            ---      NeV   2p3s
   3      NeV   2p3p           ---      NeV   2p3d
   4      NeV   2s3s           ---      NeV   2s3p
   5      NeV   2s3d           ---      NeV   2p33s
   6      NeV   2p33p          ---      NeV   2p33d
   7      NeV   2p4p           ---      NeV   2p4s
   8      NeV   2p4f           ---      NeV   2p4d
}
  if pos('ENERGIES IN UNITS',S)=0 then
    Exit;
  SamePar:=pos('DIP',S)=0;

  while pos('---',S)=0 do readln(fRCGout,S);
  if IOResult<>0 then
    Error('Error reading configurations table for transitions in RCG.'^G);
  i:=0;
  p1:=pos('---',S);
  while p1<>0 do begin
  { Read the table of configurations to find correct conf. numbers }
    inc(i);
    if i>500 then
      Error('Too many configurations.'^G);

    S1:=GetConfName(copy(S,p1-15,7));
    i1:=GetConfNo(S1);
    if (i1<1) and (S1<>'') then
      Error('Configuration '+S1+' is not found in transitions file.'^G);

    if i1>0 then begin
      i2:=i;
      ReplaceCN
    end;

    if not SamePar then begin
      S1:=GetConfName(copy(S,p1+15,7));
      i1:=GetConfNo(S1);
      if (i1<1) and (S1<>'') then
        Error('Configuration '+S1+' is not found in transitions file.'^G);

      if i1>0 then begin
        i2:=1001-i;
        ReplaceCN
      end
    end;

    readln(fRCGout,S);
    if IOResult<>0 then
      Error('Error reading configurations table for transitions in RCG.'^G);
    p1:=pos('---',S)
  end
end;

type
  BoolFunc = function(P1,P2: Pointer): boolean;
function TrLess(P1,P2: pointer): boolean;
var Pu,Pl: PLevel; El,Eu: single;
begin
  with PTrans(P1)^ do begin
    Pu:=FromLevel;
    Pl:=ToLevel
  end;
  El:=abs(Pu^.GetBestE-Pl^.GetBestE);
  with PTrans(P2)^ do begin
    Pu:=FromLevel;
    Pl:=ToLevel
  end;
  Eu:=abs(Pu^.GetBestE-Pl^.GetBestE);
  TrLess:=Eu<El
end;

function LevLess1(P1,P2: pointer): boolean;
begin
  LevLess1:=PLevel(P1)^.ET>PLevel(P2)^.ET { descending theor. energy order }
end;

function LevLess2(P1,P2: pointer): boolean;
var PL1: PLevel absolute P1;
  PL2: PLevel absolute P2;
  PC1,PC2: PConf;
  i1,i2: integer;
  S1,S2: string[3];
begin
  i1:=PAmpl(PL1^.Ampls^.At(0))^.ConfNo;
  i2:=PAmpl(PL2^.Ampls^.At(0))^.ConfNo;
  PC1:=CT[i1];
  PC2:=CT[i2];

  if i1<500 then
    i1:=-i1
  else
    i1:=1001-i1;
  if i2<500 then
    i2:=-i2
  else
    i2:=1001-i2;
  str(i1,S1);
  str(i2,S2);

  if UseNames then
    LevLess2:=(PC1^.Name<PC2^.Name) or
              (PC1^.Name=PC2^.Name) and
               ((PL1^.g<PL2^.g) or
                (PL1^.g=PL2^.g) and (PL1^.GetBestE<Pl2^.GetBestE))
  else
    LevLess2:=(S1<S2) or
              (i1=i2) and
               ((PL1^.g<PL2^.g) or
                (PL1^.g=PL2^.g) and (PL1^.GetBestE<Pl2^.GetBestE))
end;

procedure Sort(var C: TCollection; LessFunc: BoolFunc);
  { Batcher's parallel sorting algorithm }
var Pn,Qn,Qn0,Rn,Dn: LongInt;
  i,j,j1: integer;
  P1,P2: pointer;
begin
  Pn:=1;
  repeat
    Pn:=Pn*2
  until Pn>=C.Count;
  Qn0:=Pn div 2;
  repeat
    Pn:=Pn div 2;
    Qn:=Qn0;
    Rn:=0;
    Dn:=Pn;
    repeat
      j:=0;
      with C do for i:=0 to Count-dn-1 do if i and pn=rn then begin
        i1:=i+dn;
        P1:=Items^[i];
        P2:=Items^[i1];
        if LessFunc(P2,P1) then begin
          Items^[i]:=P2;
          Items^[i1]:=P1;
        end
      end;
      if qn<>pn then begin
        dn:=qn-pn;
        qn:=qn div 2;
        rn:=pn;
        j:=1
      end
    until j=0
  until Pn=1
end;

procedure FillTrans;
var i: integer;
procedure CheckLev(P: PLevel);
var MS: boolean;
procedure CheckTrans(Tr: PTrans); far;
var PR: PRange; i,j: integer;
begin
  lambda:=1e5/abs(Tr^.FromLevel^.GetBestE-Tr^.ToLevel^.GetBestE);
  if MS then
    Tr^.gA:=Tr^.gA*Fms;
  if Tr^.gA>=1e-10 then begin
    for i:=1 to Ranges.Count do begin
      PR:=Ranges.At(i-1);
      if (lambda>=PR^.FromWL) and (lambda<=Pr^.ToWL) then begin
        Tr^.gA:=Tr^.gA*PR^.gAScale/PR^.gAmax;
        if Tr^.gA>=Imin then begin
          if Trans.Count>=MaxCollectionSize then
            Error('Too many transitions.'^G);
          Trans.Insert(Tr);
          if MS then
            inc(NTms)
        end else
          Tr^.gA:=0.0;
        Exit
      end
    end;
    Tr^.gA:=0.0
  end
end;
begin { CheckLev }
  MS:=(P^.Lft>=Lms) and (abs(Fms-1)>1e-2);
  P^.AA^.ForEach(@CheckTrans);
end;
begin { FillTrans }
  for i:=0 to RCGLevels.Count-1 do
    CheckLev(RCGLevels.At(i))
end;

function LName(PL: PLevel): string;
var S1: string[11]; j: integer; PAm: PAmpl;
const
  Parity: array[boolean] of char =('_','~');
begin
  PAm:=PL^.Ampls^.At(0);
  PAm:=PAm^.TermRef;
  if UseNames then begin
    i:=PAm^.ConfNo;
    S1:=CT[i]^.Name + Parity[i>500] + PAm^.TermName
  end else begin
    i:=PAm^.ConfNo;
    if i>500 then
      j:=1001-i
    else
      j:=i;
    Str(j,S1);
    while length(S1)<2 do S1:=S1+' ';
    S1:=S1 + Parity[i>500] + '(' + copy(PAm^.TermName,1,2) + ') '
                                + copy(PAm^.TermName,3,3)
  end;
  if (length(PAm^.TermName) < 5) and not TabDelimitedOutput then
     S1 := S1 + ' ';
  LName:=S1
end;

procedure CheckCloseFileError(s: str15);
begin
  if IOResult<>0 then
    Error('Error closing '+s+' file. May be, insufficient space on disk.'^G)
end;

procedure CheckFileWriteError;
begin
  if IOResult<>0 then
    Error('File write error. May be, insufficient space on disk.'^G)
end;

const
  Qstr: array[boolean] of string[1]=(' ','?');
  Estr: array[boolean] of string[1]=(' ','*');
  Mstr: array[boolean] of string[1]=(' ','M');
procedure WriteTr(var f: text; Tr: PTrans);
var PAm,PAm1: PAmpl; LftStr: string[10]; Lft,WL: double; 
  gAstr, CFstr, gAstr1, Es1, Es2, Js1, Js2,gfstr,WLstr: str15;
begin
  with Tr^ do begin
    PAm:=ToLevel^.Ampls^.At(0);
    PAm:=PAm^.TermRef;
    PAm1:=FromLevel^.Ampls^.At(0);
    PAm1:=PAm1^.TermRef;
    Lft:=FromLevel^.Lft;
    if (Lft>=1e-6) or (Lft<0.5e-11) then
      Str(Lft*1e9:10,LftStr)
    else
      Str(Lft*1e9:10:2,LftStr);

		Str(cf:5:2,CFstr);

    WL:=1e3*abs(FromLevel^.GetBestE-ToLevel^.GetBestE);

    if (WL>MinAirWN) and (WL<MaxAirWN) then
      WL:=Lair(1e8/WL)
    else
      WL:=1e8/WL;
    if gA>=0.1 then begin
      Str(gA0:12,gAstr);
      gAstr:=copy(gAstr,1,8)+copy(gAstr,11,2);

      if TabDelimitedOutput then begin
        Trim(LftStr);
        Trim(CFstr);
        Trim(gAstr);
        Str(WL:11:3,WLstr);
        Trim(WLstr);
        Str(ToLevel^.GetBestE:11:4,Es1);
        Trim(Es1);
        Str(FromLevel^.GetBestE:11:4,Es2);
        Trim(Es2);
        Str((ToLevel^.g-1)/2:4:1,Js1);
        Trim(Js1);
        Str((FromLevel^.g-1)/2:4:1,Js2);
        Trim(Js2);
        Str(gA:7:0,gAstr1);
        Trim(gAstr1);
        Str(gf:7:3,gfstr);
        Trim(gfstr);

        writeln(f,Es1,#9,
               LName(ToLevel),#9,
               Js1,#9,
               Es2,#9,
               LName(FromLevel),#9,
               Js2,#9,WLstr,#9,gAstr1,#9,gAstr,#9,gfstr,#9,
               Qstr[Q],Estr[ToLevel^.Fixed],Estr[FromLevel^.Fixed],
               Mstr[(FromLevel^.Lft>=Lms) and (Fms>1)],#9,LftStr,#9,CFstr);
      end else
        writeln(f,ToLevel^.GetBestE:11:4,' ',
               LName(ToLevel),
               (ToLevel^.g-1)/2:4:1,' -',
               FromLevel^.GetBestE:11:4,' ',
               LName(FromLevel),
               (FromLevel^.g-1)/2:4:1,WL:11:3,gA:7:0,gAstr,gf:7:3,' ',
               Qstr[Q],Estr[ToLevel^.Fixed],Estr[FromLevel^.Fixed],
               Mstr[(FromLevel^.Lft>=Lms) and (Fms>1)],LftStr,CFstr);
    end;
  end;
  CheckFileWriteError;
end;

procedure WriteLev(var f: text; i: integer; PL: PLevel; Meta: boolean);
var dE: Real; j: byte; Es, dEs, Jstr: string[12];
begin
  if not PL^.Fixed then
    dE:=dEth
  else
    dE:=dEexp;
  if not TabDelimitedOutput then begin
    write(f,i:3);
    if not Meta then
      write(f,PL^.ET*1000:12:3,dE:10:3);
    write(f,abs(PL^.E)*1000:12:3,' ',Estr[PL^.Fixed]);
    if not Meta then
      write(f,(abs(PL^.E)-PL^.ET)*1000:12:3);
    write(f,{'  J=',}(PL^.g-1)/2:5:1,' /',LName(PL),'/   ',
          Mstr[(PL^.Lft>=Lms) and (Fms>1)]);
    {if PL^.Lft>Lms then}
    if (PL^.Lft>=1e-6) or (PL^.Lft<0.5e-11) then
      write(f,PL^.Lft*1e9:10)
    else
      write(f,PL^.Lft*1e9:10:2);
    if PL^.gList<>nil then
      for j:=1 to length(PL^.gList^) do begin
        write(f,(ord(PL^.gList^[j])-1)/2:5:1);
        if j<length(PL^.gList^) then
          write(f,',')
      end;
  end else begin
    write(f,i);
    Str(PL^.ET*1000:12:3,Es);
    Trim(Es);
    Str(dE:10:3,dEs);
    Trim(dEs);
    if not Meta then
      write(f,#9,Es,#9,dEs);
    Str(abs(PL^.E)*1000:12:3,Es);
    Trim(Es);
    write(f,#9,Es,#9,Estr[PL^.Fixed]);
    Str((PL^.g-1)/2:5:1,Jstr);
    Trim(Jstr);
    Str((abs(PL^.E)-PL^.ET)*1000:12:3,Es);
    Trim(Es);
    if not Meta then
      write(f,#9,Es);
    write(f,#9,Jstr,#9,LName(PL){,#9,
          Mstr[(PL^.Lft>=Lms) and (Fms>1)]});
    
    if (PL^.Lft>=1e-6) or (PL^.Lft<0.5e-11) then
      Str(PL^.Lft*1e9:10,Es)
    else
      Str(PL^.Lft*1e9:10:2,Es);
    Trim(Es);
    write(f,#9,Es);
    if PL^.gList<>nil then begin
      Es := '';
      for j:=1 to length(PL^.gList^) do begin
        Str((ord(PL^.gList^[j])-1)/2:5:1,Jstr);
        Es := Es + Jstr;
        if j<length(PL^.gList^) then
          Es := Es + ',';
      end;
      write(f,#9,Es);
    end;
  end;
  writeln(f);
  CheckFileWriteError;
end;

procedure WriteMetastables;
var i: integer;
procedure WriteLine(Tr: PTrans); far;
begin
  if Tr^.FromLevel^.Lft>Lms then
    WriteTr(fm,Tr)
end;
procedure WriteL(PL: PLevel); far;
begin
  inc(i);
  if PL^.Lft>Lms then
    WriteLev(fm,i,PL,true)
end;
begin
  i:=0;
  assign(fm,'meta');
  rewrite(fm);
  if IOResult<>0 then
    Error('Error creating META file.'^G);
  writeln(fm,'Metastable levels: LFT>=',Lms:9);
  writeln(fm,'Min. rad. lifetime =',Lmin*1e9:10, ' ns');
  writeln(fm,Nms,' metastable levels found.');
  if IOResult<>0 then
    Error('Error writing to META file.'^G);
  RCGLevels.ForEach(@WriteL);
  writeln(fm);
  writeln(fm,'Transitions down from metasteble levels: gA is multiplied by ',round(Fms));
  writeln(fm,'                                                       Lambda (A) gA*',round(Fms),' (s-1)');
  if IOResult<>0 then
    Error('Error writing to META file.'^G);
  Trans.ForEach(@WriteLine);
  close(fm);
  CheckCloseFileError('META');
end;

procedure WriteENLEV;
var i: integer;
procedure WriteL(PL: PLevel); far;
begin
  inc(i);
  PL^.N:=i;
  WriteLev(f4,i,PL,false)
end;
begin { WriteENLEV }
  i:=0;
  RCGLevels.ForEach(@WriteL);
  close(f4);
  CheckCloseFileError('ENLEV');
end;

procedure WriteSLFile;
procedure WriteL(PL: PLevel); far;
var Estr,Jstr: string[12];
begin
  if TabDelimitedOutput then begin
    Str(PL^.GetBestE*1000:12:3,Estr);
    Str((PL^.g-1)/2:4:1,Jstr);
    Trim(Estr);
    Trim(Jstr);
    WriteLn(f5,PL^.N,#9,Estr,#9,Jstr,#9,LName(PL));
  end else
    WriteLn(f5,PL^.N:5,PL^.GetBestE*1000:12:3,(PL^.g-1)/2:4:1,' ',LName(PL));
  if IOResult<>0 then
    Error('Error writing CONF-J-E-sorted levels file.'^G)
end;
begin { WriteSLFile }
  RCGLevels.ForEach(@WriteL);
  close(f5);
  CheckCloseFileError('sorted levels');
end;

procedure WriteTrans;
procedure WriteLine(Tr: PTrans); far;
begin
  WriteTr(f2,Tr)
end;
begin
  if not TabDelimitedOutput then begin
    writeln(f2,'          1st level                   2nd level         Lambda(A)  Int      gA    log(gf)   LFT_upp (ns) CF');
    if IOResult<>0 then
      Error('Error writing to reduced linelist file.'^G);
  end;
  Trans.ForEach(@WriteLine);
  close(f2);
  CheckCloseFileError('reduced trans');
end;

procedure WriteTRANS_DAT;
var i: integer;
procedure WriteL(PL: PLevel);
procedure WriteLower(Tr: PTrans); far;
var LInt,LNo: integer; PLin: PLine;
  WnObs,WnCalc,dWn: double;
begin  { WriteLower}
  if Tr^.gA<1e-10 then
    Exit;
  with Tr^ do begin
    LInt:=0;
    LNo:=0;
    WnObs:=0.0;
    dWn:=0.0;
    WnCalc:=(abs(FromLevel^.E)-abs(ToLevel^.E))*1e3;
    if Paint and FromLevel^.Fixed and ToLevel^.Fixed then begin
      PLin:=Tr^.ObsLine;
      if PLin<>nil then begin
        LInt:=PLin^.Int;
        LNo:=Lines.IndexOf(PLin)+1;
        WnObs:=PLin^.Wn;
        dWn:=WnObs-WnCalc
      end
    end;
    writeln(f3,'+',ToLevel^.N:4,gA:5:0,abs(ToLevel^.E)*1000:12:3,
      ' ',Estr[ToLevel^.Fixed],WnCalc:13:3,LInt:6,WnObs:12:3,dWn:11:3,LNo:6);
    if IOResult<>0 then
      Error('Error writing TRANS.DAT file.'^G)
  end
end;
begin { WriteL }
  writeln(f3,'$',i:4,'    J=',(PL^.g-1)/2:4:1,'  ',LName(PL),
          abs(PL^.E)*1000:13:3,' ',Estr[PL^.Fixed],' ');
  if IOResult<>0 then
    Error('Error writing TRANS.DAT file.'^G);
  PL^.AA^.ForEach(@WriteLower)
end;

procedure AssignDraft(PT: PTrans);
var WnUnc,WnCalc,dW0: double; LN,j: integer;
procedure Est(PL: PLine); far;
var dW: double; PN,PN1: PNumSetLine; j: integer;
begin
  inc(LN);
  PN:=nil;
  for j:=1 to Numset.Count do begin
    PN1:=Numset.At(j-1);
    if LN>=PN1^.Start then
      PN:=PN1
  end;
  if PN=nil then
    Error('Error in line numbers in NUMSET.DAT.'^G);

  WnUnc:=1e-8*PN^.dWl*sqr(PL^.Wn);
  dW:=abs(PL^.Wn-WnCalc);
  if (dW<=WnUnc) and ((PT^.ObsLine=nil) or (dW<dW0)) then begin
    PT^.ObsLine:=PL;
    dW0:=dW
  end;
end;
begin { AssignDraft }
  if (not PT^.FromLevel^.Fixed) or (not PT^.ToLevel^.Fixed) then
    Exit;
  WnCalc:=(abs(PT^.FromLevel^.E)-abs(PT^.ToLevel^.E))*1e3;
  LN:=0;
  Lines.ForEach(@Est)
end;

procedure AssignProof(PT: PTrans);
var j,j1: integer; PL,PL1: PLine; PT1: PTrans;
procedure Correct(PT1: PTrans); far;
begin
  if (PT1=PT) or (PT1^.ObsLine=nil) then
    Exit;
  if PT1^.ObsLine=PT^.ObsLine then
    if PT1^.gA/PT^.gA<AThresh then
      PT1^.ObsLine:=nil
    else
      if PT^.gA/PT1^.gA<AThresh then
        PT^.ObsLine:=nil
end;
begin { AssignProof }
  if PT^.ObsLine=nil then
    Exit;
  Trans.ForEach(@Correct)
end;
begin  { WriteTRANS_DAT }
  for i:=1 to Trans.Count do
    AssignDraft(Trans.At(i-1));
  for i:=1 to Trans.Count do
    AssignProof(Trans.At(i-1));

  for i:=1 to RCGLevels.Count do
    WriteL(RCGLevels.At(i-1));
  close(f3);
  CheckCloseFileError('TRANS.DAT');
end;

function TrLess1(P1,P2: pointer): boolean; far;
begin
  TrLess1:=LevLess1(PTrans(P1)^.ToLevel,PTrans(P2)^.ToLevel)
end;

function AmpLess(P1,P2: pointer): boolean; far;
begin
  { descending amplitude order }
  AmpLess:=abs(PAmpl(P1)^.Ampl)>abs(PAmpl(P2)^.Ampl)
end;

procedure SortTrans;
var i: integer;
procedure Sort1(PL: PLevel); far;
begin
  Sort(PL^.AA^,TrLess1)
end;
begin { SortTrans }
  RCGLevels.ForEach(@Sort1)
end;

procedure DeleteUClevels;
var PL,PL1: PLevel; i: integer;
function Unconnected(PL: PLevel): boolean;
function Found(Tr: PTrans): boolean; far;
begin
  Found:=Tr^.ToLevel=PL
end;
function Visible(Tr: PTrans): boolean; far;
begin
  Visible:=Tr^.gA>1e-10
end;
begin
  Unconnected:=((PL^.AA^.Count=0) or (PL^.AA^.FirstThat(@Visible)=nil))
                and (Trans.FirstThat(@Found)=nil)
end;
begin { DeleteUClevels }
  ND:=0;
  repeat
    PL:=nil;
    for i:=0 to RCGLevels.Count-1 do begin
      PL1:=RCGLevels.At(i);
      if Unconnected(PL1) then begin
        PL:=PL1;
        i:=RCGLevels.Count-1
      end
    end;
    if PL<>nil then begin
      RCGLevels.Delete(PL);
      Dispose(PL,Done);
      inc(ND)
    end else
      exit
  until false
end;

function TestTr(PL1,PL: PLevel): PTrans;
function TrExists(Tr: PTrans): boolean; far;
begin
  TrExists:=Tr^.ToLevel=PL
end;
begin
  TestTr:=PL1^.AA^.FirstThat(@TrExists);
end;

procedure InsertEigenRefs(Pu: PLevel; Ipar,Ij: integer);
var IL1,IL2,i: integer;
procedure EigenRefs(PAm: PAmpl); far;
var
  Nc,j: integer; PAm1: PAmpl; OK: boolean;
  Nc_str, Ipar_str: string[3];
begin
  Nc:=PAm^.ConfNo;
  OK:=false;
  for j:=IL1 to IL2 do begin
    PAm1:=LSBasis[Ipar].At(j-1);
    if (PAm1<>nil) and (PAm1^.ConfNo=Nc)
    and (PAm1^.TermName=PAm^.TermName) then begin
      PAm^.TermRef:=PAm1;
      OK:=true;
      break;
    end
  end;
  if not OK then begin
    Str(Nc, Nc_str);
    Str(Ipar, Ipar_str);
    Error('RCG term reference not found for RCE levels.' + #13#10 +
      'Ipar = ' + Ipar_str + '; Conf = ' + CT[Nc]^.Name +
      '; Term = ' + PAm^.TermName + '.' + ^G)
  end
end;
begin  { InsertEigenRefs }
  IL1:=0;
  for i:=1 to Ij-1 do
    inc(IL1,NLevJ[i,IPar]);
  IL2:=IL1+NLevJ[Ij,Ipar];
  inc(IL1);
  Pu^.Ampls^.ForEach(@EigenRefs)
end;

procedure ReadRCELevels;
var A: integer;
begin
  Stop:=eof(fRCEOUT);
  i:=0;
  i1:=0;
  Nlev:=0;
  Stop1:=false;
  El:=-1e20;
  Header:=true;
  Ipar:=1;
  NJVal[Ipar]:=0;
  gprev:=0;

  while not Stop do begin
    readln(fRCEOUT,S);
    if (S<>'') and (S[7]='.') and (S[21]='.') then begin
      Header:=false;
      if Stop1 then begin
           { first parity ended, second began }
        i1:=0;
        Ipar:=2;
        NJVal[Ipar]:=0;
        gprev:=0;
      end;
      Stop1:=false;
      inc(i);
      if MaxAvail<SizeOf(TLevel)+2*(SizeOf(TCollection)+1000) then
        Error('Insufficient memory.'^G);

      Pu:=New(PLevel,Init);
      if Pu=nil then
        Error('Low memory.'^G);

      S1:=copy(S,1,13);
      Trim(S1);
      val(S1,Pu^.E,Res);
      if Res<>0 then
        Error('Levels source file format error.'^G);

      S1:=copy(S,15,10);
      Trim(S1);
      val(S1,Pu^.ET,Res);

      if Res<>0 then
        Error('Levels source file format error.'^G);

      EL:=Pu^.ET;
      if S[14]<>'*' then
        Pu^.Fixed:=true; { An experimental datum is available for this level }

      p1:=0;
      j1:=12;
      repeat
        inc(p1);
        inc(j1,24);
        S2:=copy(S,j1,6);
        if (S2<>'') then begin
          S1:=GetConfName(S2);
          Nc:=GetConfNo(S1);
          if Nc<1 then begin
            inc(i1);
            if i1>500 then
              Error('Too many configurations.'^G);
            if MaxAvail<SizeOf(TConf)+1000 then
              Error('Insufficient memory.'^G);
            PC:=New(PConf,Init);
            if PC=nil then
              Error('Low memory.'^G);
            PC^.Name:=S1;
            if Ipar=2 then
              Nc:=1001-i1
            else
              Nc:=i1;
            CT[Nc]:=PC;
          end;

          S1:=Copy(S,j1-4,3);
          Trim(S1);
          val(S1,A,Res);
          if Res<>0 then
            Error('RCEOUT file format error for amplitude value '+S1+^G);

            { RCEOUT contains only amplitudes >= 1% }
          if MaxAvail<SizeOf(TAmpl)+1000 then
            Error('Insufficient memory.'^G);
          PAm:=New(PAmpl,Init);
          if PAm=nil then
            Error('Low memory.'^G);
          PAm^.ConfNo:=Nc;
          PAm^.TermName:=Copy(S,j1+8,2)+Copy(S,j1+12,3);
          Trim(PAm^.TermName);
          PAm^.TermRef:=nil;
          PAm^.Ampl:=A*100;
          if Pu^.Ampls^.Count>=MaxCollectionSize then
            CollectionSizeError('vector amplitudes');
          Pu^.Ampls^.Insert(PAm)
        end
      until (S2='') or (p1>=MaxAmplComps);

      S1:=copy(S,27,4);
      Pu^.g:=GetG1(S1,Res);
      if Pu^.g<>gprev then begin
        inc(NJVal[Ipar]);
        NlevJ[NJVal[Ipar],Ipar]:=0
      end;
      gprev:=Pu^.g;
      inc(NlevJ[NJVal[Ipar],Ipar]);
      if Res<>0 then
        Error('RCEOUT file format error for J value '+S1+^G);

      if RCELevels.Count >= MaxCollectionSize then
        CollectionSizeError('RCE levels list');

      RCELevels.Insert(Pu)

    end { if (S[7]='.') ... }
    else begin
      if not Header then
        Stop1:=true;
      Header:=true
    end;
    if Stop1 and (Ipar=2) or eof(fRCEOUT) then
      Stop:=true
  end;
  Nlev:=i;
  close(fRCEOUT);
end;

procedure ReadRCGEigenValues;
var i,j,Ij,IL,IL1,IL2,NJv,Ipar: integer;
  Sc,St: string;
begin
  for Ipar:=1 to 2 do begin
    for NJv:=1 to NJVal[Ipar] do begin
      { Read eigenvalues from RCG output file }
      repeat
        readln(fRCGout,S);
        if IOResult<>0 then
          Error('Error reading RCG output file.'^G);
      until eof(fRCGout) or (pos('EIGENVALUES',S)<>0);
      if pos('EIGENVALUES',S)=0 then
        Error('EIGENVALUES section missing in RCG output file.'^G);
      for Ir:=1 to (NLevJ[NJv,Ipar]-1) div 11 +1 do begin
        readln(fRCGout,S);
        readln(fRCGout,S);
        if IOResult<>0 then
          Error('Error reading RCG output file.'^G);

        IL1:=1+11*(Ir-1);
        IL2:=IL1+10;
        if IL2>NLevJ[NJv,Ipar] then
          IL2:=NLevJ[NJv,Ipar];
        for IL:=IL1 to IL2 do begin
          if MaxAvail<SizeOf(TLevel)+2*(SizeOf(TCollection)+1000) then
            Error('Insufficient memory.'^G);
          Pu:=New(PLevel,Init);
          if Pu=nil then
            Error('Low memory.'^G);
          S2:=copy(S,29+(IL-IL1)*9,9);
          Trim(S2);
          val(S2,Pu^.ET,Res);
          if (Res<>0) then
            Error('RCG output file format error in EIGENVALUES section: '+
                   S2+^G);
          Pu^.E:=Pu^.ET;

          j:=0;
          if IPar=2 then
            for i:=1 to NJVal[1] do
              inc(j,NLevJ[i,1]);
          for i:=1 to Njv-1 do
            inc(j,NLevJ[i,Ipar]);
          inc(j,IL);
          Pl:=PLevel(RCELevels.At(j-1));
          if Pl=nil then
            Error('Level search failed in RCG.'^G);
          Pu^.g:=Pl^.g;

          if (IL>1) and (Pu^.g<>gprev) then
            Error('Mismatch of number of levels in J-group between RCG and RCE.'^G);
          gprev:=Pu^.g;

          if RCGLevels.Count >= MaxCollectionSize then
            CollectionSizeError('RCG levels list');

          RCGLevels.Insert(Pu)
        end
      end;

      { Read eigenvectors from RCG OUTG11 file }
      repeat
        readln(fRCGout,S);
        if IOResult<>0 then
          Error('Error reading RCG output file.'^G);
      until eof(fRCGout) or (pos('EIGENVECTORS',S)<>0);
      if pos('EIGENVECTORS',S)=0 then
        Error('EIGENVECTORS section missing in RCG output file.'^G);
            { Initialize the basis states }
      for i:=1 to NLevJ[NJv,Ipar] do begin
        if MaxAvail<SizeOf(TAmpl)+1000 then
          Error('Insufficient memory.'^G);
        PAm:=New(PAmpl,Init);
        if PAm=nil then
          Error('Low memory.'^G);

        if LSBasis[IPar].Count >= MaxCollectionSize then
          CollectionSizeError('basis set states list');

        LSBasis[IPar].Insert(PAm);
      end;

      readln(fRCGout,S);
      if IOResult<>0 then
        Error('Error reading RCG output file.'^G);
      for Ir:=1 to (NLevJ[NJv,Ipar]-1) div 11 +1 do begin
        readln(fRCGout,S);
        readln(fRCGout,Sc);
        readln(fRCGout,St);
        if IOResult<>0 then
          Error('Error reading RCG output file.'^G);

        IL1:=1+11*(Ir-1);
        IL2:=IL1+10;
        if IL2>NLevJ[NJv,Ipar] then
          IL2:=NLevJ[NJv,Ipar];

        for Ij:=1 to NLevJ[NJv,Ipar] do begin
          readln(fRCGout,S);
          if IOResult<>0 then
            Error('Error reading RCG output file.'^G);

          i:=LSBasis[Ipar].Count-NLevJ[NJv,Ipar]+Ij;
          PAm1:=LSBasis[Ipar].At(i-1);
          if PAm1=nil then
            Error('NIL pointer in LS basis.'^G);
          if Ir=1 then begin
            S2:=copy(S,5,9);
            CName:=GetConfName(S2);
            Nc:=GetConfNo(CName);
            if (Nc<0) or ((Nc<=500) and (Ipar=2))
            or ((Nc>500) and (Ipar=1)) then
              Error('Configuration name '+S2
                    +' not identified in RCG output file.'^G);
            PAm1^.ConfNo:=Nc;
            PAm1^.TermRef:=NIL;

            PAm1^.TermName:=copy(S,18,2)+copy(S,22,3);
            Trim(PAm1^.TermName);
            PAm1^.Ampl:=10000;

            {writeln(Ipar, ':',Nc,':', PAm1^.TermName,':');}
          end;

          for IL:=IL1 to IL2 do begin
            S2:=copy(S,30+(IL-IL1)*9,8);
            Trim(S2);
            val(S2,Aj,Res);
            if Res<>0 then
              Error('RCG output file format error for component amplitude: '
                    +S2+^G);
            if abs(Aj) >= 0.005 then begin
              if MaxAvail<SizeOf(TAmpl)+1000 then
                Error('Insufficient memory.'^G);
              PAm:=New(PAmpl,Init);
              if PAm=nil then
                Error('Low memory.'^G);
              j:=0;
              if IPar=2 then
                for i:=1 to NJVal[1] do
                  inc(j,NLevJ[i,1]);
              for i:=1 to Njv-1 do
                inc(j,NLevJ[i,Ipar]);
              inc(j,IL);
              Pu:=RCGLevels.At(j-1);
              if Pu=nil then
                Error('Error in procedure reading EIGENVECTORS.'^G);

              PAm^.Ampl:=round(Aj*10000);
              PAm^.TermRef:=PAm1;
              S2:=copy(Sc,31+(IL-IL1)*9,9);
              CName:=GetConfName(S2);
              Nc:=GetConfNo(CName);
              if (Nc<0) or ((Nc<=500) and (Ipar=2))
              or ((Nc>500) and (Ipar=1)) then
                Error('Configuration name '+S2
                      +' not identified in RCG output file.'^G);
              PAm^.ConfNo:=Nc;
              PAm^.TermName:=copy(St,32+(IL-IL1)*9,2)+
                             copy(St,36+(IL-IL1)*9,3);
              Trim(PAm^.TermName);
              if Pu^.Ampls^.Count>=MaxCollectionSize then
                CollectionSizeError('vector amplitudes');
              Pu^.Ampls^.Insert(PAm)
            end
          end
        end
      end;
    end
  end;
end;

procedure IdentifyRCGLevels;
var i,j,j1,NJv,Ir,Ipar,IRev: integer;
begin
  IL:=0;
  for Ipar:=1 to 2 do begin
    for Njv:=1 to NJval[Ipar] do begin
      IL1:=IL+1;
      NL:=NLevJ[Njv,Ipar];
      IL2:=IL+NL;
      for Ir:=1 to NL do begin
        inc(IL);
        Pu:=RCGLevels.At(IL-1);
        MinDist:=100000000;
        for Irev:=1 to 2 do begin
          for i:=IL1 to IL2 do begin
            Pl:=RCELevels.At(i-1);
            Dist:=0;
            jm:=MaxAmplComps;
            if jm>NL then
              jm:=NL;
            Nfound:=0;
            for j:=1 to jm do begin
              if j>Pu^.Ampls^.Count then
                break;
              PAm:=Pu^.Ampls^.At(j-1);
              A1:=PAm^.Ampl;
              if Irev=2 then
                A1:=-PAm^.Ampl; { Check for inverted sign }
              A2:=0;
              for j1:=1 to jm do begin
                if j1>Pl^.Ampls^.Count then
                  break;
                PAm1:=Pl^.Ampls^.At(j1-1);
                if PAm1^.TermRef=PAm^.TermRef then begin
                  inc(Nfound);
                  A2:=PAm1^.Ampl;
                  break
                end
              end;
              Dist:=Dist+sqr(LongInt(A1-A2));
            end;
            if Nfound<jm then begin { Add components of Pl missing in Pu }
              for j:=1 to jm do begin
                if j>Pl^.Ampls^.Count then
                  break;
                PAm1:=Pl^.Ampls^.At(j-1);
                for j1:=1 to jm do begin
                  if j1>Pu^.Ampls^.Count then
                    break;
                  PAm:=Pu^.Ampls^.At(j1-1);
                  if PAm^.TermRef=PAm1^.TermRef then
                    break
                end;
                if PAm^.TermRef<>PAm1^.TermRef then
                  Dist:=Dist+sqr(LongInt(PAm1^.Ampl))
              end
            end;
            if Dist<MinDist then begin
              MinDist:=Dist;
              BestMatch:=i
            end
          end
        end;
        PL:=RCELevels.At(BestMatch-1);
        if Pu^.N<>0 then
          Error('Failed to identify RCG eigenvectors with RCE eigenvectors'^G);
        {with PL^ do
          if (E>=0) and (ET>=0) or (E<=0) and (ET<=0) then}
        if PL^.Fixed then
          Pu^.E:=PL^.E;
        Pu^.Fixed:=PL^.Fixed;
        PL^.ET:=Pu^.ET;
        { Pu^.ET:=PL^.ET; }
        Pu^.N:=BestMatch
        {Pu^.ConfNo:=PL^.ConfNo
         Pu^.Term:=PL^.Term}
      end
    end
  end;
end;

procedure ReadTransitions;
var gf: real;
procedure Add_gA(PT: PTrans); far;
begin
  Pu^.Lft:=Pu^.Lft+PT^.gA
end;
begin
  Start:=false;
  n:=0;
  Stop:=eof(fRCGout);
  Pu:=nil;

  while not eof(fRCGout) do begin
    readln(fRCGout,S);
    if IOResult<>0 then
      Error('Error reading transitions file.'^G);
    ReadCTable;
    if not Start then begin
      p1:=pos('LAMBDA',S);
      Start:=p1<>0
    end else begin
      Stop:=eof(fRCGout);
      if pos('* * *',S)<>0 then begin
        Pu:=nil;
        if not GetLevel(S,Pu) then begin
{         Pu:=nil;
          close(f2);
          halt(1)}

        end
      end;
      p1:=pos('E+',S);
      if (p1<>0) and (pos('SUMF',S)=0) and (pos('GAMIN',S)=0)
      and (pos('GM*',S)=0) and (pos('GA=',S)=0) then begin
{
        if Pu=nil then begin
          close(f2);
          Error('Undefined initial level in transitions file.'^G)
        end;
}
        S1:=copy(S,p1-5,9);
        val(S1,gA,Res);

        if Res=0 then begin
          S1:=copy(S,p1-12,6);
          val(S1,gf,Res);
        end;
        if Res=0 then begin
          S1:=Copy(S,p1+5,5);
          Trim(S1);
          val(S1,Br,Res)
        end;
        if (Res<>0) and (N=0) then
          Error('Transitions source file format error.');
        if not GetLevel(S,Pl) then begin
          Pl:=nil;
{
          close(f2);
          halt
}
        end;
        if (Pu<>nil) and (Pl<>nil) then begin
          if {(Pu^.ET>Pl^.ET) or not (KPrint in [3,7])} true then begin
            if Pu^.ET>Pl^.ET then
              PL1:=Pu
            else
              PL1:=PL; { PL1 always signifies level with higher energy }

            if PL1<>Pu then
              PL:=Pu; { PL is lower energy level }

            Tr:=TestTr(PL1,PL);
            E1:=PL1^.GetBestE;
            E2:=PL^.GetBestE;
            dE:=abs(E1-E2);
            if (gA>=1e-10) and (Tr=nil) and (dE>1e-10)
              and (E1<=Emax) and (E2<=Emax)
            then begin
              if MaxAvail<sizeof(TTrans)+1000 then begin
                close(f2);
                Error('Insufficient memory.'^G)
              end;
              lambda:=1e5/dE;
              Tr:=New(PTrans,Init);
              if Tr=nil then
                Error('Low memory.'^G);
              Tr^.FromLevel:=PL1;
              Tr^.ToLevel:=PL;
              Tr^.gA:=gA;
              Tr^.gA0:=gA;
              Tr^.gf:=gf;
							Tr^.cf:=abs(Br);
              Tr^.Q:=abs(Br)<0.1;
              Tr^.ObsLine:=nil;
              if PL1^.AA^.Count>=MaxCollectionSize then
                CollectionSizeError('transitions from a single level');
              PL1^.AA^.Insert(Tr);
              {Trans.Insert(Tr);}
              inc(n)
            end
          end
        end
      end
    end
  end;
  close(fRCGout);
  Nms:=0;
  for i:=1 to RCGLevels.Count do begin
    Pu:=RCGLevels.At(i-1);
    Pu^.Lft:=0.0;
    Pu^.AA^.ForEach(@Add_gA);
    Pu^.Lft:=Pu^.Lft/Pu^.g;
    if Pu^.Lft>1e-20 then
      Pu^.Lft:=1/Pu^.Lft
    else
      Pu^.Lft:=1e20;
    if Pu^.Lft>Lms then
      inc(Nms);
    if Pu^.Lft<Lmin then
      Lmin:=Pu^.Lft
  end
end;

procedure MergeTerms;
var i,j,j0,k,Nmerged: integer; {PAm,PAm1: PAmpl;} Tr,Tr1: PTrans;
  PL,Pu,PT: PLevel;
function SameLower(PT1: PTrans): boolean; far;
begin
  SameLower:=PT1^.ToLevel=PT
end;

begin { MergeTerms }
  write('Merging terms with splitting <= ',
         round(MinTermSplit*1e3),' cm-1 ... ');
  Nmerged:=0;
  i:=0;
  while i<RCGLevels.Count-1 do begin
    Pu:=RCGLevels.At(i);
    PAm:=Pu^.Ampls^.At(0);
    j:=0;
    j0:=-1;
    while j<RCGLevels.Count do begin
      if j=i then begin
        inc(j);
        continue
      end;
      PL:=RCGLevels.At(j);
      PAm1:=PL^.Ampls^.At(0);
      if {(PL^.g<>Pu^.g) and}
        (PAm^.TermRef^.ConfNo=PAm1^.TermRef^.ConfNo)
      and (PAm^.TermName=PAm1^.TermName)
      and ((abs(Pu^.E-PL^.E)<=MinTermSplit) or
         (abs(Pu^.ET-PL^.ET)<=MinTermSplit)
        and (abs(PAm^.Ampl)>=5500)
        and (abs(PAm1^.Ampl)>=5500)) then
      begin
        while PL^.AA^.Count>0 do begin  { Transfer transitions from PL to Pu }
          Tr:=PL^.AA^.At(0);
          PT:=Tr^.ToLevel;
          Tr1:=Pu^.AA^.FirstThat(@SameLower);
          if Tr1<>nil then begin        { Add to Tr1 if found }
            if (Tr1^.Q) and (not Tr^.Q) and { Handle 'questionable' flag }
              (Tr^.gA>=Tr1^.gA) then begin
              Tr1^.Q:=false;
							Tr1^.cf:=Tr^.cf; { Accept the max. cf value }
						end;
            Tr1^.gA:=Tr1^.gA+Tr^.gA;
            Tr1^.gA0:=Tr1^.gA0+Tr^.gA0;
            Tr1^.gf:=Tr1^.gf+Tr^.gf;
            {if Nmerged>=98 then
              write('   *1: ');}
            Dispose(Tr,Done)
            {if Nmerged>=98 then
              writeln('passed.');}
          end else begin       { insert in Pu^.AA^ transitions list }
            Tr^.FromLevel:=Pu;
            {if Nmerged>=98 then begin
              if Pu^.AA=nil then
                writeln('AA = nil pointer.');
              writeln('   *2: MaxAvail=',MaxAvail, '; AA^.Count=',Pu^.AA^.Count,'. ');
            end;}
            if PL1^.AA^.Count>=MaxCollectionSize then
              CollectionSizeError('transitions from a single level');
            Pu^.AA^.Insert(Tr);
            {if Nmerged>=98 then
              writeln('passed.');}
          end;
          {if Nmerged>=98 then
            write('   *3: ');}
          PL^.AA^.AtDelete(0);
          {if Nmerged>=98 then
            writeln('passed.');}
        end;
        for k:=0 to RCGLevels.Count-1 do
          if (k<>i) and (k<>j) then begin
            PL1:=RCGLevels.At(k);
            PT:=PL;
            Tr:=PL1^.AA^.FirstThat(@SameLower);
            if Tr<>nil then begin
              PT:=Pu;
              Tr1:=PL1^.AA^.FirstThat(@SameLower);
              if Tr1<>nil then begin
                if (Tr1^.Q) and (not Tr^.Q) and
                  (Tr^.gA>=Tr1^.gA) then begin
                  Tr1^.Q:=false;
		  Tr1^.cf:= Tr^.cf; { Accept the max. cf value }
		end;
                Tr1^.gA:=Tr1^.gA+Tr^.gA;
                Tr1^.gA0:=Tr1^.gA0+Tr^.gA0;
                Tr1^.gf:=Tr1^.gf+Tr^.gf;
                {if Nmerged>=98 then
                  write('   *4 ');}
                PL1^.AA^.Delete(Tr);
                Dispose(Tr,Done);
                {if Nmerged>=98 then
                  writeln('passed.');}
              end else
                TR^.ToLevel:=Pu
            end
          end;
        Pu^.ET:=(Pu^.ET*Pu^.g+PL^.ET*PL^.g)/(Pu^.g+PL^.g);
        if Pu^.Fixed then
          if PL^.Fixed then
            Pu^.E:=(Pu^.E*Pu^.g+PL^.E*PL^.g)/(Pu^.g+PL^.g)
          else
            Pu^.E:=Pu^.E
        else
          if PL^.Fixed then
            Pu^.E:=PL^.E
          else
            Pu^.E:=-Pu^.ET;
        if Pu^.gList=nil then begin
          if MaxAvail<SizeOf(Pu^.gList)+1000 then
            Error('Not enough memory.'^G);
          GetMem(Pu^.gList,SizeOf(Pu^.gList^));
          fillchar(Pu^.gList^,SizeOf(Pu^.gList^),#0);
          inc(Pu^.gList^[0]);
          Pu^.gList^[length(Pu^.gList^)]:=chr(Pu^.g);
        end;
        Pu^.g:=Pu^.g+PL^.g;
        if PL^.gList=nil then begin
          inc(Pu^.gList^[0]);
          if ord(Pu^.gList^[0])>MaxJval then
            Error('Too many J values.'^G);
          Pu^.gList^[length(Pu^.gList^)]:=chr(PL^.g);
        end else
          for k:=1 to length(PL^.gList^) do begin
            inc(Pu^.gList^[0]);
            if ord(Pu^.gList^[0])>MaxJval then
              Error('Too many J values.'^G);
            Pu^.gList^[length(Pu^.gList^)]:=PL^.gList^[k]
          end;
        {if Nmerged>=98 then
          write('   *5 ');}
        RCGLevels.Delete(PL);
        Dispose(PL,Done);
        {if Nmerged>=98 then
          writeln('passed.');}
        inc(Nmerged);
        {writeln(Nmerged, i:5, j:5);}
        if j<i then
          dec(i);
        if j0<0 then
          j0:=j
      end else
        inc(j)
    end;
    if (j0<i) and (j0>=0) then
      i:=j0
    else
      inc(i)
  end;
  if Nmerged=1 then
    writeln('done. 1 level is deleted.')
  else
    writeln('done. ',Nmerged,' levels are deleted.')
end;

procedure ReadLines;
var PL: PLine;
begin
  if not Lines.Init(100,10) then
    Error('Low memory.'^G);
  assign(fDLV,'DLV.DAT');
  reset(fDLV);
  if IOResult<>0 then
    Error('Error opening DLV.DAT file.'^G);
  while not eof(fDLV) do begin
    readln(fDLV,S);
    if IOResult<>0 then
      Error('Error reading DLV.DAT file.'^G);
    if MaxAvail<SizeOf(TLine)+1000 then
      Error('Insufficient memory.'^G);
    PL:=New(PLine,Init);
    if PL=nil then
      Error('Low memory.'^G);
    S2:=copy(S,1,5);
    Trim(S2);
    val(S2,PL^.Int,Res);
    if Res<>0 then
      Error('Format error in line intensity in DLV.DAT: '+S2+^G);
    S2:=copy(S,6,14);
    Trim(S2);
    val(S2,PL^.Wn,Res);
    if Res<>0 then
      Error('Format error in line wavenumber in DLV.DAT: '+S2+^G);
    {PL^.Comments:=copy(S,20,length(S));}
    if Lines.Count>=MaxCollectionSize then
      Error('Too many lines.'^G);
    Lines.Insert(PL);
  end;
  close(fDLV);
end;

procedure ReadNumset;
var
  PN: PNumsetLine;
begin
  if not NumSet.Init(100,10) then
    Error('Low memory.'^G);
  assign(fNUMSET,'NUMSET.DAT');
  reset(fNUMSET);
  if IOResult<>0 then
    Error('Error opening NUMSET.DAT file.'^G);
  while not eof(fNUMSET) do begin
    readln(fNUMSET,S);
    if IOResult<>0 then
      Error('Error reading NUMSET.DAT file.'^G);
    if S[1]<>'$' then begin
      if MaxAvail<SizeOf(TNumsetLine)+1000 then
        Error('Insufficient memory.'^G);
      PN:=New(PNumsetLine,Init);
      if PN=nil then
        Error('Low memory.'^G);
      S2:=copy(S,2,5);
      Trim(S2);
      val(S2,PN^.Start,Res);
      if Res<>0 then
        Error('Format error in group start number in NUMSET.DAT.'^G);
      if (PN^.Start<=0) or (PN^.Start>Lines.Count) then
        Error('Group start number out of range in NUMSET.DAT.'^G);
      S2:=copy(S,24,5);
      Trim(S2);
      val(S2,PN^.dWl,Res);
      if Res<>0 then
        Error('Format error in wavelength uncertainty in NUMSET.DAT.'^G);
      Numset.Insert(PN)
    end
  end;
  close(fNUMSET)
end;

procedure Find_gAmax;
procedure Test(PT: PTrans); far;
var PR: PRange; i: integer; lambda: real;
begin
  lambda:=1e5/abs(PT^.FromLevel^.GetBestE-PT^.ToLevel^.GetBestE);
  if PT^.ToLevel=PGround then
    ResTrans.Insert(PT);
  for i:=1 to Ranges.Count do begin
    PR:=Ranges.At(i-1);
    if (lambda>=PR^.FromWL) and (lambda<=Pr^.ToWL) then begin
      if not PR^.Absorption then begin
        if PR^.Teff>1e-10 then
          PT^.gA:=PT^.gA*exp(-PT^.FromLevel^.ET/PR^.Teff)/lambda;
        if (PT^.gA>PR^.gAmax) then
          PR^.gAmax:=PT^.gA
      end;
      break
    end
  end
end;
procedure Test1(PT: PTrans); far;
var A: real;
begin
  A:=PT^.gA/PT^.FromLevel^.g;
  if A>AgrMax then
    AgrMax:=A;
end;
procedure Test2(PT: PTrans); far;
var PR: PRange; i,j: integer; lambda: real; TR1: PTrans; PL,PL1: PLevel;
begin
  lambda:=1e5/abs(PT^.FromLevel^.GetBestE-PT^.ToLevel^.GetBestE);
  for i:=1 to Ranges.Count do begin
    PR:=Ranges.At(i-1);
    if (lambda>=PR^.FromWL) and (lambda<=Pr^.ToWL) then begin
      if PR^.Absorption then begin
        if PT^.ToLevel=PGround then
          PL:=PGround
        else
          PL:=nil;
        if PL=nil then for j:=1 to ResTrans.Count do begin
           TR1:=ResTrans.At(j-1);
           PL1:=TR1^.FromLevel;
           if PL1=PT^.ToLevel then begin
             PL:=PL1;
             break
           end
        end;
        if PL<>nil then begin
            { Lower level is populated and transition can be observed
              in absorption }
          {PT^.gA:=(PT^.gA*PT^.ToLevel^.g/PT^.FromLevel^.g)*
                   (TR1^.gA/TR1^.FromLevel^.g)/AgrMax}
          PT^.gA:=(PT^.gA*PL^.g/PT^.FromLevel^.g);
          if PR^.Teff>1e-10 then
            PT^.gA:=PT^.gA*exp(-PL^.ET/PR^.Teff)
        end else
          PT^.gA:=0.0;
        if PT^.gA>PR^.gAmax then
          PR^.gAmax:=PT^.gA;
      end;
      break
    end
  end
end;
var i: integer; PL: PLevel;
begin { Find_gAmax }
  for i:=1 to RCGLevels.Count do begin
       { Scan emission lines }
    PL:=RCGLevels.At(i-1);
    PL^.AA^.ForEach(@Test)
  end;
  AgrMax:=0;
  ResTrans.ForEach(@Test1);
  for i:=1 to RCGLevels.Count do begin
       { Scan absorption lines }
    PL:=RCGLevels.At(i-1);
    PL^.AA^.ForEach(@Test2)
  end
end;

begin { Main program }
   { Initialize global collections }
  writeln('Memory available: ',MaxAvail);
  if not RCGLevels.Init(500,10) then
    Error('Low memory.'^G);
  if not RCELevels.Init(500,10) then
    Error('Low memory.'^G);
  if not Trans.Init(500,10) then
    Error('Low memory.'^G);
  if not LSBasis[1].Init(100,10) then
    Error('Low memory.'^G);
  if not LSBasis[2].Init(100,10) then
    Error('Low memory.'^G);
  if not Ranges.Init(10,1) then
    Error('Low memory.'^G);
  if not ResTrans.Init(100,10) then
    Error('Low memory.'^G);
     { Read input/output parameters from file CONV_OUT.PAR }
  ReadParams;
     { Initialize variables }
  fillchar(CT,sizeof(CT),0);
  Nms:=0;
  Lmin:=1e20;

  { Read levels list RCEOUT created by RCE }
  write('Reading RCE levels ... ');
  ReadRCELevels;
  writeln('done.');


  write('Reading RCG eigenvalues and eigenvectors ... ');
  ReadRCGEigenValues;
  writeln('done.');

     { Sort RCG eigenvectors }
  IL:=0;
  for Ipar:=1 to 2 do
    for Ij:=1 to NJVal[Ipar] do begin
      for i:=1 to NLevJ[Ij,Ipar] do begin
        inc(IL);
        Pu:=RCGLevels.At(IL-1);
        if Pu=nil then
          Error('NIL pointer in RCGLevels.'^G);
        Pu^.N:=0;
        Sort(Pu^.Ampls^,AmpLess)
      end;
          { Check for duplicates }
      IL1:=IL-NLevJ[Ij,Ipar]+1;
      IL2:=IL;
      for i:=IL1 to IL2-1 do begin
        Pu:=RCGLevels.At(i-1);
        PAm:=Pu^.Ampls^.At(0);
        for j:=i+1 to IL2 do begin
          PL:=RCGLevels.At(j-1);
          PAm1:=PL^.Ampls^.At(0);
          if (abs(Pu^.ET-PL^.ET)<0.8e-3) and
          (PAm^.TermRef=PAm1^.TermRef) then
            writeln('Warning: duplicate level energy',Pu^.ET:9:3,' with J=',
                     (Pu^.g-1)/2:4:1,' in RCG.'^G);
        end
      end
    end;

      { Insert term references into RCE eigenvectors }
  IL:=0;
  for Ipar:=1 to 2 do
    for Ij:=1 to NJVal[Ipar] do
      for i:=1 to NLevJ[Ij,Ipar] do begin
        inc(IL);
        Pu:=RCELevels.At(IL-1);
        if Pu=nil then
          Error('NIL pointer in RCELevels.'^G);
          InsertEigenRefs(Pu,Ipar,Ij)
      end;

      { Identify RCG levels with RCE levels }
  write('Identifying RCG levels with RCE levels ... ');
  IdentifyRCGLevels;
  writeln('done.');
  NLms:=0;
  NTms:=0;
    { Read all transitions from RCG OUTG11 file }
  write('Reading transitions data ... ');
  ReadTransitions;
  writeln('done.');

  if MinTermSplit>0.8e-3 then
    { Merge terms with small splittings into single levels }
    MergeTerms;

  { Sort the levels by energy }
  write('Sorting levels ... ');
  Sort(RCGLevels,LevLess1);
  writeln('done.');
  PGround:=RCGLevels.At(RCGLevels.Count-1);

  Find_gAmax;

     { Write the header of the reduced linelist file }
  for i:=1 to Ranges.Count do begin
    PR:=Ranges.At(i-1);
    if not TabDelimitedOutput then
      writeln(f2,'Range',i:2,':',PR^.FromWL:10:3,'-',PR^.ToWL:10:3,
            ': Max. gA =',PR^.gAmax:10, ' sec-1, Scale =',PR^.gAScale:5:0,
            ': Teff =',PR^.Teff*1000:12:0,' cm-1 : ',Absorp[PR^.Absorption]);
    if IOResult<>0 then
      Error('Error writing to reduced linelist file.'^G)
  end;

    {  Reduce the linelist: }
    {  write only strong transitions to reduced file }
  NTms:=0;
  write('Reducing transitions file ... ');
  FillTrans;
  Sort(Trans,TrLess);
  WriteTrans;
  if (Fms>1) and (NTms<>0) then
    WriteMetastables;
  writeln('done.');

  write('Sorting transitions ... ');
  SortTrans;
  writeln('done.');

  if RemoveUnconnected then begin
    ND:=0; i:=0;
    write('Removing unconnected levels ... ');
    DeleteUClevels;
    if ND=1 then
      writeln('done. 1 level is removed.')
    else
      writeln('done. ',ND,' levels are removed.')
  end;

  assign(f4,'enlev.dat');
  rewrite(f4);
  if IOResult<>0 then
    Error('Error creating ENLEV.DAT file.'^G);

  write('Writing levels to ENLEV.DAT ... ');
  WriteENLEV;
  writeln('done.');

  if Paint then begin
    ReadLines;
    ReadNumset
  end;

  assign(f3,'trans.dat');
  rewrite(f3);
  if IOResult<>0 then
    Error('Error creating TRANS.DAT file.'^G);

  write('Writing trans. probabilities to TRANS.DAT ... ');
  WriteTRANS_DAT;
  writeln('done.');

  if MakeSLFile then begin
    write('Creating levels file sorted by CONF-J-E order ... ');
    Sort(RCGLevels,LevLess2);
    WriteSLFile
  end;
  writeln('done.');

  writeln;
  for i:=1 to 2 do
    LSBasis[i].Done;
  RCGLevels.Done;
  RCELevels.Done;
  Trans.DeleteAll;
  Trans.Done;
  Ranges.Done;
  if Paint then begin
    NumSet.Done;
    Lines.Done
  end;
  writeln('Ok.');
end.