uses crt,graph;
const
sysN=4;
at=1.0;
alf1=0.01;bet1=1.0;
alf2=1.0;bet2=0.01;
r=1.0;
RES=5000; {chislo shagov Runge-Kutta za period T}
xmax=1.5;xmin=-1.5;
ymax=1.5;ymin=-1.5;

CRE=-1.0; CIM=-0.0;
DRE=1.0; DIM=0.0;

type
matr1=array[1..sysN] of extended;

var
gx,gy,gd,gm:integer;
dyn:matr1;
jjj1,jjj:LongInt;
phi,phi1:extended;

{============================================================================}
{============================================================================}

procedure dynamics(f:matr1; alf,bet,t:extended; var f_:matr1);
begin
f_[1]:=alf*(f[1]*f[4]-f[2]*f[3]);
f_[2]:=-alf*(f[1]*f[3]+f[2]*f[4]);
f_[3]:=bet*(2*f[1]*f[2]);
f_[4]:=-bet*(f[1]*f[1]-f[2]*f[2]);
end;
{============================================================================}
{============================================================================}

Procedure RungeKutta(h,alf,bet,t:extended; var pf:matr1);
var i_:integer;
k:array[1..4,1..sysN]of extended;
f,f_:matr1;
begin
for i_:=1 to sysN do begin f[i_]:=pf[i_]; end;
dynamics(f,alf,bet,t,f_);
for i_:=1 to sysN do begin k[1,i_]:=h*f_[i_]; end;

for i_:=1 to sysN do begin f[i_]:=pf[i_]+k[1,i_]/2; end;
dynamics(f,alf,bet,t+h/2,f_);
for i_:=1 to sysN do begin k[2,i_]:=h*f_[i_]; end;

for i_:=1 to sysN do begin f[i_]:=pf[i_]+k[2,i_]/2; end;
dynamics(f,alf,bet,t+h/2,f_);
for i_:=1 to sysN do begin k[3,i_]:=h*f_[i_]; end;

for i_:=1 to sysN do begin f[i_]:=pf[i_]+k[3,i_]; end;
dynamics(f,alf,bet,t+h,f_);
for i_:=1 to sysN do begin k[4,i_]:=h*f_[i_]; end;

for i_:=1 to sysN do begin pf[i_]:=pf[i_]+(1/6)*(k[1,i_]+2*k[2,i_]+2*k[3,i_]+k[4,i_]); end;
end;
{============================================================================}
{============================================================================}
BEGIN
gd:=detect;gm:=3; InitGraph(gd,gm,'c:\FPC\2.2.0\units\i386-win32\graph');

rectangle(0,0,300,300);rectangle(300,0,600,300);

dyn[1]:=0.10; dyn[2]:=0.0; dyn[3]:=0.0; dyn[4]:=0.0;phi:=0;
FOR jjj:=1 to 5000 DO BEGIN
for jjj1:=1 to RES do RungeKutta(r/RES,alf1,bet1,jjj1*r/RES,dyn);
dyn[1]:=DRE; dyn[2]:=DIM;
dyn[3]:=dyn[3]/sqrt(sqr(dyn[3])+sqr(dyn[4])); dyn[4]:=dyn[4]/sqrt(sqr(dyn[3])+sqr(dyn[4]));
for jjj1:=1 to RES do RungeKutta(r/RES,alf2,bet2,jjj1*r/RES,dyn);
dyn[1]:=dyn[1]+CRE; dyn[2]:=dyn[2]+CIM; dyn[3]:=0.0; dyn[4]:=0.0;
    gx:=round(300*(dyn[1]-xmin)/(xmax-xmin));gy:=round(300*(dyn[2]-ymin)/(ymax-ymin));putpixel(gx,300-gy,white);
   if dyn[1]>0 then phi1:=arctan(dyn[2]/(dyn[1]+10e-10)) else phi1:=pi+arctan(dyn[2]/(dyn[1]+10e-10));if phi1<0 then phi1:=phi1+2*pi;
   gx:=round(300*phi/6.28);gy:=round(300*phi1/6.28);putpixel(300+gx,300-gy,white);   phi:=phi1;
dyn[1]:=at*dyn[1]/sqrt(sqr(dyn[1])+sqr(dyn[2])); dyn[2]:=at*dyn[2]/sqrt(sqr(dyn[1])+sqr(dyn[2]));
END;

readln;
END.
