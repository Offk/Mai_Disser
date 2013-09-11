unit Unit1;

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, StdCtrls, ExtCtrls;

type
  TForm1 = class(TForm)
    Button1: TButton;
    KDRImage: TImage;
    procedure Button1Click(Sender: TObject);
    procedure FormClose(Sender: TObject; var Action: TCloseAction);

  private
    { Private declarations }
  public
    { Public declarations }
  end;


var
  Form1: TForm1;
  KDR,Lyap: TBitmap;


const
  ReAmin=-1.5;
  ReAmax=1.5;
  ImAmin=-1.5;
  ImAmax=1.5;
  PixSizeA=300;
  PixSizeB=300;

  alf1=0.1;
  bet1=1;
  alf2=1;
  bet2=0.1;
  R1=1;
  R2=1;
  ReD=1;
  ImD=0;
  ReC=-1.5;
  ImC=0;

  sysN=4;
  eps=0.001;
  VeryBigNumber=10000;
type
  vector = array[1..sysN] of double;
  vectArr = array[1..sysN] of vector;
implementation

{$R *.dfm}

function f(x:vector; alf,bet:double):vector;
begin
  f[1]:=alf*(x[1]*x[4]-x[2]*x[3]);
  f[2]:=-alf*(x[1]*x[3]+x[2]*x[4]);
  f[3]:=bet*(2*x[1]*x[2]);
  f[4]:=-bet*(x[1]*x[1]-x[2]*x[2]);
end;

function df(x:vector; dx:vector; alf,bet:double):vector;
begin
  df[1]:=alf*(x[1]*dx[4]+dx[1]*x[4]-x[2]*dx[3]-dx[2]*x[3]);
  df[2]:=-alf*(x[1]*dx[3]+dx[1]*x[3]+x[2]*dx[4]+dx[2]*x[4]);
  df[3]:=bet*(2*x[1]*dx[2]+2*dx[1]*x[2]);
  df[4]:=-bet*(x[1]*dx[1]+dx[1]*x[1]-x[2]*dx[2]-dx[2]*x[2]);
end;

function f_runge(x:vector; R:double; alf,bet:double):vector;
var
  h:double;
  t,x1,x2,x3,x4,res:vector;
  k:integer;
begin

  h:=R/100;

  x1:=f(x,alf,bet);

  for k:=1 to sysN do
    t[k]:=x[k]+h*x1[k]/2;

  x2:=f(t,alf,bet);

  for k:=1 to sysN do
    t[k]:=x[k]+h*x2[k]/2;

  x3:=f(t,alf,bet);

  for k:=1 to sysN do
    t[k]:=x[k]+h*x3[k];

  x4:=f(t,alf,bet);

  for k:=1 to sysN do
    res[k]:=x[k]+h*(x1[k]+2*x2[k]+2*x3[k]+x4[k])/6;

  result:=res;
end;

function df_runge(x:vector; dx:vector; R:double; alf,bet:double):vector;
var
  h:double;
  t,x1,x2,x3,x4,res:vector;
  k:integer;
begin

  h:=R/100;

  x1:=df(x,dx,alf,bet);

  for k:=1 to sysN do
    t[k]:=dx[k]+h*x1[k]/2;

  x2:=df(x,t,alf,bet);

  for k:=1 to sysN do
    t[k]:=dx[k]+h*x2[k]/2;

  x3:=df(x,t,alf,bet);

  for k:=1 to sysN do
    t[k]:=dx[k]+h*x3[k];

  x4:=df(x,t,alf,bet);

  for k:=1 to sysN do
    res[k]:=dx[k]+h*(x1[k]+2*x2[k]+2*x3[k]+x4[k])/6;

  result:=res;
end;

function matrEvol(x:vector; dx:vectArr; R:double; alf,bet:double):vectArr;
var res:vectArr;
    k: integer;
begin
  for k:=1 to sysN do
    res[k]:=df_runge(x,dx[k],R,alf,bet);
  result:=res;
end;

function length(x:vector):double;
var res:double;
    k:integer;
begin
  res:=0;
  for k:=1 to sysN do
    res:=res+x[k]*x[k];
  result:=sqrt(res);
end;

function vectArr_length(x:vectarr):double;
var res:double;
    k:integer;
begin
  res:=0;
  for k:=1 to sysN do
    res:=res+length(x[k]);
  result:=sqrt(res);
end;

function max(x,y:double):double;
begin
  if (x>y)
    then
      max:=x
    else
      max:=y;
end;

function mult(x:vector; a:double):vector;
var res:vector;
    k:integer;
begin
  for k:=1 to sysN do
    res[k]:=res[k]*a;
  result:=res;
end;

function plus(x,y:vector):vector;
var res:vector;
    k:integer;
begin
  for k:=1 to sysN do
    res[k]:=x[k]+y[k];
  result:=res;
end;

function NullVector():vector;
var k:integer;
    res:vector;
begin
  for k:=1 to sysN do
    res[k]:=0;
  result:=res;
end;

function RunawayVector():vector;
var k:integer;
    res:vector;
begin
  for k:=1 to sysN do
    res[k]:=999;
  result:=res;
end;

function norm(x:vector):vector;
var res:vector;
begin
  if length(x)>0.0001
    then
      res:=mult(x,1/length(x))
    else
      begin
        res:=NullVector();
      end;
  result:=res;
end;

function scal(x,y:vector):double;
var k:integer;
    res:double;
begin
  res:=0;
  for k:=1 to sysN do
    res:=res+x[k]*y[k];
  result:=res;
end;

function ElementaryPertuberation():vectArr;
var res:vectArr;
    k: integer;
begin
  for k:=1 to sysN do
    res[k]:=NullVector();
  for k:=1 to sysN do
    res[k][k]:=1;
  result:=res;
end;

procedure ortogonalization(dx:vectArr);
var i,j:integer;
    temp:vector;
begin
  for i:=1 to sysN do
    begin
      j:=2;
      while (j<i) do
        begin
          temp:=mult(mult(dx[i],scal(dx[j],dx[i])),-1);
          dx[i]:=plus(dx[i],temp);
          inc(j);
        end;
      dx[i]:=norm(dx[i]);
    end;
end;

procedure depth_cross(R,alf,bet:double; var x:vector; var runaway:boolean);
var j:integer;
begin
  for j:=1 to 100 do
    begin
      x:=f_Runge(x,R,alf,bet);
      if (length(x)>VeryBigNumber) then
        begin
          runaway:=true;
          exit;
        end;
    end;
end;


function per_circle(x:vector):vector;
var module_a,module_b:double;
    runaway:boolean;
begin
  runaway:=false;
  x[3]:=0;
  x[4]:=0;
  module_a:=sqrt(x[1]*x[1]+x[2]*x[2]);
  x[1]:=x[1]/module_a;
  x[2]:=x[2]/module_a;

  depth_cross(R1,alf1,bet1,x,runaway);
  if runaway then
    begin
      result:=RunawayVector();
      exit;
    end;
  x[1]:=ReD;
  x[2]:=ImD;

  module_b:=sqrt(x[3]*x[3]+x[4]*x[4]);
  //if module_b<0.00001 then
  //  module_b:=0.00001;
  x[3]:=x[3]/module_b;
  x[4]:=x[4]/module_b;

  depth_cross(R2,alf2,bet2,x,runaway);
  if runaway then
    begin
      result:=RunawayVector();
      exit;
    end;
  x[1]:=x[1]+ReC;
  x[2]:=x[2]+ImC;
  //x[3]:=0;
  //x[4]:=0;
  result:=x;
end;

procedure SetPixel_PhaseProtret(x,y:double);
var i,j:integer;
begin
  i:=round(PixSizeA*(x-ReAmin)/(ReAmax-ReAmin));
  j:=round(PixSizeB-PixSizeB*(y-ImAmin)/(ImAmax-ImAmin));
  with KDR.Canvas do
  begin
    Pixels[i,j]:=ClBlack;
  end;
end;

function Iterate(x0:vector):vector;
const
  PredCicleNum=1000;
  CicleNum=10000;
var i:integer;
    x:vector;
begin
  x:=x0;
  {
  for i:=0 to PredCicleNum do
  begin
    x:=per_circle(x);
    if x[1]=999 then
      exit;
  end;
  }
  for i:=0 to CicleNum do
  begin
    x:=per_circle(x);
    if x[1]=999 then
      exit;
    SetPixel_PhaseProtret(x[3],x[2]);
        form1.KDRImage.Canvas.CopyRect(Rect(0,0,PixSizeA,PixSizeB),KDR.canvas,Rect(0,0,PixSizeA,PixSizeB));
        application.ProcessMessages;
  end;

  result:=x;

end;

procedure TForm1.Button1Click(Sender: TObject);
var x:vector;
    i:integer;
begin
  KDR:=TBitmap.Create;
  KDR.Width:=PixSizeA;
  KDR.Height:=PixSizeB;
  with KDR.Canvas do
    begin


      for i:=1 to 1 do
      begin
        x[1]:=-1.5+Random*3.0;
        x[2]:=-1.5+Random*3.0;
        x[3]:=-1.5+Random*3.0;
        x[4]:=-1.5+Random*3.0;

        x:=Iterate(x);

      end;
    end;

  DateSeparator:= '-';
  TimeSeparator:= '_';
  KDR.SaveToFile(DateTimeToStr(Now)+'.bmp');
end;

procedure TForm1.FormClose(Sender: TObject; var Action: TCloseAction);
begin
   halt(1);
end;

end.

