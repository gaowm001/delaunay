#ifndef DELAUNAY_H
#define DELAUNAY_H

#include <stdio.h>
#include <math.h>
#include <string>
#include <QDebug>
#include <omp.h>

typedef struct {
    int more=-1;
    int less=-1;
    int sort=-1;
} chaindot;

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

double calcos(int a,int b,int c,double **dot,double *f) {//求<BAC
    double abx=dot[b][0]-dot[a][0],
           aby=dot[b][1]-dot[a][1],
           acx=dot[c][0]-dot[a][0],
           acy=dot[c][1]-dot[a][1];
    double r=acos((abx*acx+aby*acy)/sqrt((abx*abx+aby*aby)*(acx*acx+acy*acy)));
    *f=abx*acy-acx*aby;
    return r;
}

int isincircle(int a,int b,int c,int d,double** dot) {//d在a,b,c的圆内
    double x0,y0;
    double x1=dot[a][0],x2=dot[b][0],x3=dot[c][0],y1=dot[a][1],y2=dot[b][1],y3=dot[c][1];
    double xx=2*((y1-y2)*(x1-x3)-(x1-x2)*(y1-y3));
    if (xx==0) {
        qDebug()<<"line";
    }
    x0=((y1-y2)*(x1*x1-x3*x3+y1*y1-y3*y3)-(y1-y3)*(x1*x1-x2*x2+y1*y1-y2*y2))/xx;
    y0=((x1-x3)*(x1*x1-x2*x2+y1*y1-y2*y2)-(x1-x2)*(x1*x1-x3*x3+y1*y1-y3*y3))/xx;
    double r=(x1-x0)*(x1-x0) + (y1-y0)*(y1-y0);
    double r1=((dot[d][0]-x0)*(dot[d][0]-x0)+(dot[d][1]-y0)*(dot[d][1]-y0));
    //r>r1,d点在圆心内
    if (r-r1>1e-7) return 1;
    else if (r-r1<1e-7) return -1;
    return 0;
}

void choicedot(SVdot vl,SVdot vr,double **dot,int *l,int *r) {
    double f,rr;
    for (int i=0;i<vl.vlen;i++) {
        if (vl.v[i]==*l) continue;
        rr=calcos(*l,*r,vl.v[i],dot,&f);
        if (f==0||rr==0) {//三点同线，选离r点近的
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
        if (f==0||rr==0) {//三点同线，选离l点近的
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

void choicecan(Scandot *rdot,int a1,int a2,SVdot vl,double **dot,int *dl,bool l) {
    int j=0;*dl=-1;
      for (int i=0;i<vl.linelen;i++) {
          if (vl.line[i]==nullptr) continue;
          if (vl.line[i][0]==a2) rdot[i].dot=vl.line[i][1];
          else if (vl.line[i][1]==a2) rdot[i].dot=vl.line[i][0];
          else continue;
        double f;
        double r=calcos(a2,a1,rdot[i].dot,dot,&f);
        if ((l&&(f<0))||((!l)&&(f>0))||(f==0)||r==0) continue;
        rdot[i].r=r;
        if (*dl==-1) {
            *dl=i;
            rdot[i].less=-1;
            rdot[i].more=-1;
            continue;
        }
        j=*dl;
        while (true) {
            f=false;
            if (r==rdot[j].r) {//角度相同，离a2近的在前
                    if ((dot[rdot[i].dot][0]-dot[a2][0])*(dot[rdot[i].dot][0]-dot[a2][0])+(dot[rdot[i].dot][1]-dot[a2][1])*(dot[rdot[i].dot][0]-dot[a2][0])
                            >(dot[rdot[j].dot][0]-dot[a2][0])*(dot[rdot[j].dot][0]-dot[a2][0])+(dot[rdot[j].dot][1]-dot[a2][1])*(dot[rdot[j].dot][0]-dot[a2][0]))
                    f=true;
                qDebug()<<"samecircle:"<<rdot[i].dot<<","<<rdot[j].dot<<","<<a1<<","<<a2;
            }
            if ((r>rdot[j].r)||(f)) {
                if (rdot[j].more==-1) {
                    rdot[i].less=j;rdot[i].more=-1;
                    rdot[j].more=i;
                    break;
                }
                j=rdot[j].more;
            } else {
                rdot[i].more=j;rdot[i].less=rdot[j].less;
                if (rdot[j].less==-1) {
                    *dl=i;
                } else {
                    rdot[rdot[j].less].more=i;
                }
                rdot[j].less=i;
                break;
            }
        }
    }    
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

SVdot conquer(SVdot *v,SVdot *vl,SVdot *vr,double** dot) {
    v->linelen=vl->linelen+vr->linelen;
    if (v->linelen>=v->maxlinelen) {
        v->maxlinelen=v->linelen*2;
        v->line=new int*[v->maxlinelen];
    }
    v->linelen=vl->linelen+vr->linelen;
    memcpy(v->line,vl->line,(vl->linelen)*sizeof(vl->line));
    memcpy(v->line+vl->linelen,vr->line,(vr->linelen)*sizeof(vr->line));
    int head,dl=-1,dr=-1;
    int l1=vl->v[vl->vlen-1];
    int r1=vr->v[0];
    choicedot(*vl,*vr,dot,&l1,&r1);
    Scandot *rdot=new Scandot[vr->linelen],*ldot=new Scandot[vl->linelen];
    while (true) {
        choicecan(rdot,l1,r1,*vr,dot,&head,false);
        int i=head;
        while (i!=-1) {//右边
            if (rdot[i].more==-1) {dr=rdot[i].dot;break;}
            int r=isincircle(l1,r1,rdot[i].dot,rdot[rdot[i].more].dot,dot);
            if (r==-1) {
                dr=rdot[i].dot;
                break;
            } else if (r==1) {
                delline(v,r1,rdot[i].dot);
                vr->line[i]=nullptr;
            }
            i=rdot[i].more;
        }
        if (head==-1) dr=-1;
        choicecan(ldot,r1,l1,*vl,dot,&head,true);
        i=head;
        while (i!=-1) {//左边
            if (ldot[i].more==-1) {dl=ldot[i].dot;break;}
            int r=isincircle(l1,r1,ldot[i].dot,ldot[ldot[i].more].dot,dot);
            if (r==-1) {
                dl=ldot[i].dot;
                break;
            } else {
                delline(v,l1,ldot[i].dot);
                vl->line[i]=nullptr;
            }
            i=ldot[i].more;
        }
        if (head==-1) dl=-1;
        insertline(v,l1,r1);
        if ((dl==-1)&&(dr==-1)) {//无候选点，合并完成
            break;
        }
        if ((dl>-1)&&(dr>-1)) {
            int r=isincircle(r1,l1,dr,dl,dot);
            if (r==1) {//dl在圆内，选择dl
               dr=-1;
            } else if (r==-1) {
               dl=-1;
            } else {
                qDebug()<<"lrcircle:"<<l1<<","<<r1<<","<<dl<<","<<dr;
            }
        }
        if (dl==-1) r1=dr;
        else l1=dl;
    }
    delete[] ldot;
    delete[] rdot;
    delete[] vl->line;
    delete[] vr->line;
    return *v;
}

SVdot divide(SVdot *vdot,double **dot) {
    if (vdot->vlen>3) {
        SVdot v[2];
        v[0].vlen=vdot->vlen/2;
        v[0].v=vdot->v;
        v[1].vlen=vdot->vlen-v[0].vlen;
        v[1].v=vdot->v+v[0].vlen;
        #pragma omp parallel for
        for (int i=0;i<2;i++) {
            v[i]=divide(&v[i],dot);
        }
        return conquer(vdot,&v[0],&v[1],dot);
    } else {
        if (vdot->vlen==2) {            
            vdot->linelen=1;
            vdot->maxlinelen=vdot->linelen*2;
            vdot->line=new int*[vdot->maxlinelen];
            vdot->line[0]=new int[2];
            vdot->line[0][0]=vdot->v[0];
            vdot->line[0][1]=vdot->v[1];
        } else {//处理三点同线
            double f,r=calcos(vdot->v[0],vdot->v[1],vdot->v[2],dot,&f);
            if (f==0||r==0) {
//                qDebug()<<"same line";
                vdot->linelen=2;
                vdot->maxlinelen=vdot->linelen*2;
                vdot->line=new int*[vdot->maxlinelen];
                vdot->line[0]=new int[2];
                vdot->line[0][0]=vdot->v[0];
                vdot->line[0][1]=vdot->v[1];
                vdot->line[1]=new int[2];
                vdot->line[1][0]=vdot->v[1];
                vdot->line[1][1]=vdot->v[2];
            } else {
                vdot->linelen=3;
                vdot->maxlinelen=vdot->linelen*2;
                vdot->line=new int*[vdot->maxlinelen];
                vdot->line[0]=new int[2];
                vdot->line[0][0]=vdot->v[0];
                vdot->line[0][1]=vdot->v[1];
                vdot->line[1]=new int[2];
                vdot->line[1][0]=vdot->v[0];
                vdot->line[1][1]=vdot->v[2];
                vdot->line[2]=new int[2];
                vdot->line[2][0]=vdot->v[1];
                vdot->line[2][1]=vdot->v[2];
            }
        }
        return SVdot(*vdot);
    }
}

SVdot calDelaunay(double **dot) {
    int lens=_msize(dot)/sizeof(dot[0]);
    chaindot chain[lens];
    int dm=0,dl=0,dr=0,sum=lens;
    float dmm=0;
    for (int i=1;i<lens;i++) {
        int j=dm;
        int f=dcompare(dot[i],dot[j]);
        if (f==0) {
//            qDebug()<<"same position!"<<QString::number(i)<<","<<QString::number(j);
            sum--;
            continue;
        }
        if (f==1) {
            while (true) {
                if (chain[j].more==-1) {
                    chain[i].less=j;chain[i].more=-1;
                    chain[j].more=i;
                    dr=i;
                    break;
                }
                f=dcompare(dot[i],dot[chain[j].more]);
                if (f==0) {
//                    qDebug()<<"same position!"<<QString::number(i)<<","<<QString::number(j);
                    sum--;
                    break;
                }
                if (f<0) {
                    chain[i].less=j;chain[i].more=chain[j].more;
                    chain[chain[j].more].less=i;chain[j].more=i;
                    break;
                }
                j=chain[j].more;
            }
            dmm+=0.5;
        } else {
            while (true) {
                if (chain[j].less==-1) {
                    chain[i].more=j;chain[i].less=-1;
                    chain[j].less=i;
                    dl=i;
                    break;
                }
                f=dcompare(dot[i],dot[chain[j].less]);
                if (f==0) {
//                    qDebug()<<"same position!"<<QString::number(i)<<","<<QString::number(j);
                    sum--;
                    break;
                }
                if (f>0) {
                    chain[i].more=j;chain[i].less=chain[j].less;
                    chain[chain[j].less].more=i;chain[j].less=i;
                    break;
                }
                j=chain[j].less;
            }
            dmm-=0.5;
        }
        if (dmm==1) {dm=chain[dm].more;dmm=0;};
        if (dmm==-1) {dm=chain[dm].less;dmm=0;};
    }
    SVdot vdot;
    vdot.v=new int[sum];
    vdot.vlen=sum;
    vdot.linelen=0;
    vdot.maxlinelen=vdot.vlen*2;
    int i=dl,j=0;
    while (i!=-1) {
//        qDebug()<<i<<"x,y:"<<dot[i][0]<<","<<dot[i][1];
        vdot.v[j]=i;
        chain[i].sort=i;
        i=chain[i].more;
        j++;
    }
//    qDebug()<<"------------------------";
    return divide(&vdot,dot);
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
                qDebug()<<"line1("<<QString::number(a)<<","<<QString::number(b)<<"),line2("<<QString::number(v.line[j][0])<<","<<QString::number(v.line[j][1])<<")";
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
            if (r==0||f==0) {
                if (dot[c][0]>dot[a][0]&&dot[c][0]<dot[b][0]) {
                    qDebug()<<"line1("<<QString::number(a)<<","<<QString::number(b)<<"),dot("<<QString::number(c)<<")";
                    return false;
                }
            } else if (f>0) {
                if (x1==-1) x1=c;
                else if (isincircle(a,b,c,x1,dot)>=0) continue;
                else x1=c;
                if (x2==-1) continue;
                if (isincircle(a,b,c,x2,dot)>0) {
                    qDebug()<<"line1("<<QString::number(a)<<","<<QString::number(b)<<"),dot("<<QString::number(c)<<","<<QString::number(x2)<<")";
                    return false;
                }
            } else {
                if (x2==-1) x2=c;
                else if (isincircle(a,b,c,x2,dot)>=0) continue;
                else x2=c;
                if (x1==-1) continue;
                if (isincircle(a,b,c,x1,dot)>0) {
                    qDebug()<<"line2("<<QString::number(a)<<","<<QString::number(b)<<"),dot("<<QString::number(c)<<","<<QString::number(x1)<<")";
                    return false;
                }
            }
        }
    }
    return true;
}

#endif // DELAUNAY_H
