//============================================================================//
//                                                                            //
// delaunay                                                                   //
//                                                                            //
// A 2D Delaunay Triangulator                                                 //
//                                                                            //
// Author: Water Fox                                                          //
//                                                                            //
//                                                                            //
// Version 1.0.0                                                              //
// Apr 31, 2024                                                               //
//                                                                            //
//                                                                            //
// delaunay is freely available through the website:                          //
//       http://github.com/gaowm001/delaunay.                                 //
//   It may be copied, modified, and redistributed for non-commercial use.    //
//   Please consult the file LICENSE for the detailed copyright notices.      //
//                                                                            //
//============================================================================//

#ifndef DELAUNAY_H
#define DELAUNAY_H

#define iszero(x) fabs(x)<1e-12

#include <stdio.h>
#include <math.h>
#include <string>

typedef struct {
    int more=-1;
    int less=-1;
    int dot=-1;
    double r=0;
} Scandot;

typedef struct {
    int *v;//顶点
    int **line;//边
    int vlen;
    int maxlinelen=0;
    int linelen=0;
} SVdot;

int dcompare(double p1[],double p2[]) {
    if (p1[0]>p2[0]) return 1;
    if (p1[0]<p2[0]) return -1;
    if (p1[1]>p2[1]) return 1;
    if (p1[1]<p2[1]) return -1;
    return 0;
}

void quickSort(double** data, int start, int end,int* sort)  //并行快排
{
    if (start < end) {
        int pos=start,tend=end,temp = sort[pos];   //以第一个元素为基准
        while (pos < tend) {
            while ((pos < tend) && (dcompare(data[sort[tend]] , data[temp])>=0))
                tend--;   //找到第一个比基准小的数
            sort[pos] = sort[tend];
            while ((pos < tend) && (dcompare(data[sort[pos]],data[temp])<=0))
                pos++;    //找到第一个比基准大的数
            sort[tend] = sort[pos];
        }
        sort[pos] = temp;   //以基准作为分界线
        quickSort(data, start, pos - 1,sort);
        quickSort(data, pos + 1, end,sort);
    }
}

double calcos(int a,int b,int c,double **dot,double *f) {//求<BAC
    double abx=dot[b][0]-dot[a][0];
    double aby=dot[b][1]-dot[a][1];
    double acx=dot[c][0]-dot[a][0];
    double acy=dot[c][1]-dot[a][1];
    *f=abx*acy-acx*aby;
    double r=acos((abx*acx+aby*acy)/sqrt((abx*abx+aby*aby)*(acx*acx+acy*acy)));
    return r;
}

int isincircle(int a,int b,int c,int d,double** dot) {//d在a,b,c的圆内
    double x0,y0;
    double x1=dot[a][0];
    double x2=dot[b][0];
    double x3=dot[c][0];
    double y1=dot[a][1];
    double y2=dot[b][1];
    double y3=dot[c][1];
    double xx=2*((y1-y2)*(x1-x3)-(x1-x2)*(y1-y3));
    x0=((y1-y2)*(x1*x1-x3*x3+y1*y1-y3*y3)-(y1-y3)*(x1*x1-x2*x2+y1*y1-y2*y2))/xx;
    y0=((x1-x3)*(x1*x1-x2*x2+y1*y1-y2*y2)-(x1-x2)*(x1*x1-x3*x3+y1*y1-y3*y3))/xx;
    double ra=sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0)),rb=sqrt((x2-x0)*(x2-x0) + (y2-y0)*(y2-y0)),rc=sqrt((x3-x0)*(x3-x0) + (y3-y0)*(y3-y0)),rmax,rmin;
    if (ra>rb) {rmax=ra;rmin=rb;}
    else {rmax=rb;rmin=ra;}
    if (rc>rmax) rmax=rc;
    if (rc<rmin) rmin=rc;
    double r1=sqrt((dot[d][0]-x0)*(dot[d][0]-x0)+(dot[d][1]-y0)*(dot[d][1]-y0));
    //r>r1,d点在圆心内
    if (iszero(r1-rmin)||iszero(r1-rmax)) return 0;
    if (r1>rmax) return -1;
    else if (r1<rmin) return 1;
    return 0;
}

void insertline(SVdot *v,int a1,int a2) {
    if (v->linelen>=v->maxlinelen) {
        v->maxlinelen=v->linelen*2;
        int** temp=new int*[v->maxlinelen];
        memcpy(temp,v->line,v->linelen*sizeof(v->line));
        delete[] v->line;
        v->line=temp;
    }
    v->line[v->linelen]=new int[2];
    v->line[v->linelen][0]=a1;
    v->line[v->linelen][1]=a2;
    v->linelen++;
}

void delline(SVdot *v,int a1,int a2) {
    for (int i=0;i<v->linelen;i++) {
        if (((v->line[i][0]==a1)&&(v->line[i][1]==a2))||((v->line[i][0]==a2)&&(v->line[i][1]==a1))) {
            int** temp=new int*[v->maxlinelen];
            memcpy(temp,v->line,i*sizeof(v->line));
            memcpy(temp+i,v->line+i+1,sizeof(v->line)*(v->linelen-i-1));
            delete[] v->line[i];
            delete[] v->line;
            v->line=temp;
            v->linelen--;
            break;
        }
    }
}

void choicedot(SVdot vl,SVdot vr,double **dot,int *l,int *r) {
    double f,rr;
    for (int i=0;i<vl.vlen;i++) {
        if (vl.v[i]==*l) continue;
        rr=calcos(*l,*r,vl.v[i],dot,&f);
        if (iszero(f)) {//三点同线，选离r点近的
            if ((dot[vl.v[i]][0]-dot[*r][0])*(dot[vl.v[i]][0]-dot[*r][0])+(dot[vl.v[i]][1]-dot[*r][1])*(dot[vl.v[i]][1]-dot[*r][1])
                    <(dot[*l][0]-dot[*r][0])*(dot[*l][0]-dot[*r][0])+(dot[*l][1]-dot[*r][1])*(dot[*l][1]-dot[*r][1])) {
               *l=vl.v[i];
               choicedot(vl,vr,dot,l,r);
               return;
           }
        }
        if (f<0) {
            *l=vl.v[i];
            choicedot(vl,vr,dot,l,r);
            return ;
        }
    }
    for (int i=0;i<vr.vlen;i++) {
        if (vr.v[i]==*r) continue;
        rr=calcos(*l,*r,vr.v[i],dot,&f);
        if (iszero(f)) {//三点同线，选离l点近的
            if ((dot[vr.v[i]][0]-dot[*l][0])*(dot[vr.v[i]][0]-dot[*l][0])+(dot[vr.v[i]][1]-dot[*l][1])*(dot[vr.v[i]][1]-dot[*l][1])
                    <(dot[*l][0]-dot[*r][0])*(dot[*l][0]-dot[*r][0])+(dot[*l][1]-dot[*r][1])*(dot[*l][1]-dot[*r][1])) {
               *r=vr.v[i];
               choicedot(vl,vr,dot,l,r);
               return;
           }
        }
        if (f<0) {
            *r=vr.v[i];
            choicedot(vl,vr,dot,l,r);
            return ;
        }
    }
    return;
}

int choicecanl(int a1,int a2,SVdot *vl,double **dot) {
    Scandot *ldot=new Scandot[vl->linelen];
    int j=0,dl=-1;
    for (int i=0;i<vl->linelen;i++) {
        if (vl->line[i][0]==a2) ldot[i].dot=vl->line[i][1];
        else if (vl->line[i][1]==a2) ldot[i].dot=vl->line[i][0];
        else continue;
        double f;
        double r=calcos(a2,a1,ldot[i].dot,dot,&f);
        if (iszero(f)||(f<0)) continue;
        ldot[i].r=r;
        if (dl==-1) {
            dl=i;
            ldot[i].less=-1;
            ldot[i].more=-1;
            continue;
        }
        j=dl;
        while (true) {
            f=0;
            if (iszero(r-ldot[j].r)) {//角度相同，离a2近的在前
                    if ((dot[ldot[i].dot][0]-dot[a2][0])*(dot[ldot[i].dot][0]-dot[a2][0])+(dot[ldot[i].dot][1]-dot[a2][1])*(dot[ldot[i].dot][0]-dot[a2][0])
                            >(dot[ldot[j].dot][0]-dot[a2][0])*(dot[ldot[j].dot][0]-dot[a2][0])+(dot[ldot[j].dot][1]-dot[a2][1])*(dot[ldot[j].dot][0]-dot[a2][0]))
                    f=1;
            }
            if ((r>ldot[j].r)||(f==1)) {
                if (ldot[j].more==-1) {
                    ldot[i].less=j;ldot[i].more=-1;
                    ldot[j].more=i;
                    break;
                }
                j=ldot[j].more;
            } else {
                ldot[i].more=j;ldot[i].less=ldot[j].less;
                if (ldot[j].less==-1) {
                    dl=i;
                } else {
                    ldot[ldot[j].less].more=i;
                }
                ldot[j].less=i;
                break;
            }
        }
    }
    j=dl;
    while (j!=-1) {//左边
        if (ldot[j].more==-1) {dl=ldot[j].dot;break;}
        int r=isincircle(a2,a1,ldot[j].dot,ldot[ldot[j].more].dot,dot);
        if (r==-1) {
            dl=ldot[j].dot;
            break;
        } else {
            delline(vl,a2,ldot[j].dot);
        }
        j=ldot[j].more;
    }
    delete[] ldot;
    return dl;
}

int choicecanr(int a1,int a2,SVdot *vr,double **dot) {
    Scandot *rdot=new Scandot[vr->linelen];
    int j=0,dr=-1;
    for (int i=0;i<vr->linelen;i++) {
        if (vr->line[i][0]==a2) rdot[i].dot=vr->line[i][1];
        else if (vr->line[i][1]==a2) rdot[i].dot=vr->line[i][0];
        else continue;
        double f;
        double r=calcos(a2,a1,rdot[i].dot,dot,&f);
        if (iszero(f)||(f>0)) continue;
        rdot[i].r=r;
        if (dr==-1) {
            dr=i;
            rdot[i].less=-1;
            rdot[i].more=-1;
            continue;
        }
        j=dr;
        while (true) {
            f=0;
            if (iszero(r-rdot[j].r)) {//角度相同，离a2近的在前
                    if ((dot[rdot[i].dot][0]-dot[a2][0])*(dot[rdot[i].dot][0]-dot[a2][0])+(dot[rdot[i].dot][1]-dot[a2][1])*(dot[rdot[i].dot][0]-dot[a2][0])
                            >(dot[rdot[j].dot][0]-dot[a2][0])*(dot[rdot[j].dot][0]-dot[a2][0])+(dot[rdot[j].dot][1]-dot[a2][1])*(dot[rdot[j].dot][0]-dot[a2][0]))
                    f=1;
            }
            if ((r>rdot[j].r)||(f==1)) {
                if (rdot[j].more==-1) {
                    rdot[i].less=j;rdot[i].more=-1;
                    rdot[j].more=i;
                    break;
                }
                j=rdot[j].more;
            } else {
                rdot[i].more=j;rdot[i].less=rdot[j].less;
                if (rdot[j].less==-1) {
                    dr=i;
                } else {
                    rdot[rdot[j].less].more=i;
                }
                rdot[j].less=i;
                break;
            }
        }
    }
    j=dr;
    while (j!=-1) {//右边
        if (rdot[j].more==-1) {dr=rdot[j].dot;break;}
        int r=isincircle(a1,a2,rdot[j].dot,rdot[rdot[j].more].dot,dot);
        if (r==-1) {
            dr=rdot[j].dot;
            break;
        } else if (r==1) {
            delline(vr,a2,rdot[j].dot);
        }
        j=rdot[j].more;
    }
    delete[] rdot;
    return dr;
}

SVdot divide(SVdot vdot,double **dot) {
    if (vdot.vlen>3) {
        SVdot vl,vr;
        vl.vlen=vdot.vlen/2;
        vl.v=vdot.v;
        vr.vlen=vdot.vlen-vl.vlen;
        vr.v=vdot.v+vl.vlen;
        vl=divide(vl,dot);
        vr=divide(vr,dot);
        SVdot v;
        v.vlen=vl.vlen+vr.vlen;
        v.v=vl.v;
        v.maxlinelen=(vl.linelen+vr.linelen);
        v.line=new int*[v.maxlinelen];
        int dl=-1,dr=-1;
        int l1=vl.v[vl.vlen-1];
        int r1=vr.v[0];
        choicedot(vl,vr,dot,&l1,&r1);
        while (true) {
            dl=choicecanl(r1,l1,&vl,dot);
            dr=choicecanr(l1,r1,&vr,dot);
            insertline(&v,l1,r1);
            if ((dl==-1)&&(dr==-1)) {//无候选点，合并完成
                break;
            }
            if ((dl>-1)&&(dr>-1)) {
                int r=isincircle(r1,l1,dl,dr,dot);
                if (r==1) {//dr在圆内，选择dr
                   dl=-1;
                } else if (r==-1) {
                   dr=-1;
                }
            }
            if (dl==-1) r1=dr;
            else l1=dl;
        }
        if (v.linelen+vl.linelen+vr.linelen>v.maxlinelen) {
            v.maxlinelen=(v.linelen+vl.linelen+vr.linelen);
            int **temp=new int*[v.maxlinelen];
            memcpy(temp,v.line,(v.linelen)*sizeof(v.line));
            delete[] v.line;
            v.line=temp;
        }
        memcpy(v.line+v.linelen,vl.line,(vl.linelen)*sizeof(vl.line));
        memcpy(v.line+v.linelen+vl.linelen,vr.line,(vr.linelen)*sizeof(vr.line));
        v.linelen=v.linelen+vl.linelen+vr.linelen;
        delete[] vl.line;
        delete[] vr.line;
        return v;
    } else {
        if (vdot.vlen==2) {
            vdot.linelen=1;
            vdot.maxlinelen=vdot.linelen*2;
            vdot.line=new int*[vdot.maxlinelen];
            vdot.line[0]=new int[2];
            vdot.line[0][0]=vdot.v[0];
            vdot.line[0][1]=vdot.v[1];
        } else {//处理三点同线
            double f,r=calcos(vdot.v[0],vdot.v[1],vdot.v[2],dot,&f);
            if (iszero(f)) {
                vdot.linelen=2;
                vdot.maxlinelen=vdot.linelen*2;
                vdot.line=new int*[vdot.maxlinelen];
                vdot.line[0]=new int[2];
                vdot.line[0][0]=vdot.v[0];
                vdot.line[0][1]=vdot.v[1];
                vdot.line[1]=new int[2];
                vdot.line[1][0]=vdot.v[1];
                vdot.line[1][1]=vdot.v[2];
            } else {
                vdot.linelen=3;
                vdot.maxlinelen=vdot.linelen*2;
                vdot.line=new int*[vdot.maxlinelen];
                vdot.line[0]=new int[2];
                vdot.line[0][0]=vdot.v[0];
                vdot.line[0][1]=vdot.v[1];
                vdot.line[1]=new int[2];
                vdot.line[1][0]=vdot.v[0];
                vdot.line[1][1]=vdot.v[2];
                vdot.line[2]=new int[2];
                vdot.line[2][0]=vdot.v[1];
                vdot.line[2][1]=vdot.v[2];
            }
        }
        return vdot;
    }
}

SVdot calDelaunay(double **dot) {
    int lens=_msize(dot)/sizeof(dot[0]);
    int* sort=new int[lens];
    for (int i=0;i<lens;i++) { sort[i]=i;}
    quickSort(dot,0,lens-1,sort);
    int s=0;
    int *temp=new int[lens];
    for (int i=0;i<lens-1;i++) {
        if (dcompare(dot[sort[i]],dot[sort[i+1]])==0) {s++;continue;}
        temp[i-s]=sort[i];
    }
    temp[lens-s-1]=sort[lens-1];
    SVdot vdot;
    vdot.v=new int[lens-s];
    memcpy(vdot.v,temp,(lens-s)*sizeof(temp[0]));
    vdot.vlen=lens-s;
    delete[] sort;
    delete[] temp;
//    qDebug()<<"------------------------";
    return divide(vdot,dot);
}

void delvdot(SVdot vdot) {
    delete[] vdot.v;
    for (int i=0;i<vdot.linelen;i++) delete[] vdot.line[i];
    delete[] vdot.line;
}

bool checkDelaunay(SVdot v,double** dot) {
    for (int i=0;i<v.linelen-2;i++) {
        int a=v.line[i][0],b=v.line[i][1];
        for (int j=i+1;j<v.linelen-1;j++) {
            if (v.line[j][0]==a||v.line[j][1]==a||v.line[j][0]==b||v.line[j][1]==b) continue;
            double f,
            r1=sin(calcos(v.line[j][0],a,v.line[j][1],dot,&f)),
            r2=sin(calcos(v.line[j][0],b,v.line[j][1],dot,&f)),
            r3=sin(calcos(a,v.line[j][1],b,dot,&f)),
            r4=sin(calcos(a,v.line[j][0],b,dot,&f));
            if (r1*r2<=0&&r3*r4<=0) {
                return false;
            }
        }
    }
    return true;
}

bool checkDelaunay1(SVdot v,double ** dot) {
    for (int i=0;i<v.linelen;i++) {
        int a=v.line[i][0],b=v.line[i][1];
        int x1=-1,x2=-1;
        for (int j=0;j<v.vlen;j++) {
            if (v.v[j]==a||v.v[j]==b) continue;
            int c=v.v[j];
            double f,r=calcos(a,b,c,dot,&f);
            if (iszero(f)) {
                if (dot[c][0]>dot[a][0]&&dot[c][0]<dot[b][0]) {
                    return false;
                }
            } else if (f>0) {
                if (x1==-1) x1=c;
                else if (isincircle(a,b,c,x1,dot)>=0) continue;
                else x1=c;
                if (x2==-1) continue;
                double r=isincircle(a,b,c,x2,dot);
                if (r>0) {
                    return false;
                }
            } else {
                if (x2==-1) x2=c;
                else if (isincircle(a,b,c,x2,dot)>=0) continue;
                else x2=c;
                if (x1==-1) continue;
                double r=isincircle(a,b,c,x1,dot);
                if (r>0) {
                    return false;
                }
            }
        }
    }
    return true;
}

#endif // DELAUNAY_H
