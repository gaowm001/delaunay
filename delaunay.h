#ifndef DELAUNAY_H
#define DELAUNAY_H

#include <stdio.h>
#include <math.h>
#include <string>
#include <QDebug>

typedef struct {
    int more=-1;
    int less=-1;
    int sort=-1;
} chaindot;

typedef struct {
    int more=-1;
    int less=-1;
    double r=0;
} Scandot;

typedef struct {
    int *v;//顶点
    int **line;//边
    int vlen;
    int maxlinelen;
    int linelen=0;
} SVdot;

bool dcompare(double p1[],double p2[]) {
    if (p1[0]>p2[0]) return true;
    if (p1[0]<p2[0]) return false;
    if (p1[1]>p2[1]) return true;
    if (p1[1]<p2[1]) return false;
    return NULL;
}

double calcos(int a,int b,int c,double **dot,bool *f) {//求<BAC
    *f=true;
    double abx=dot[b][0]-dot[a][0],
           aby=dot[b][1]-dot[a][1],
           acx=dot[c][0]-dot[a][0],
           acy=dot[c][1]-dot[a][1];
    double r=acos((abx*acx+aby*acy)/sqrt((abx*abx+aby*aby)*(acx*acx+acy*acy)));
    if (abx*acy-acx*aby<0) *f=false;
    return r;
}

void choicedot(SVdot vl,SVdot vr,double **dot,int *l,int *r) {
    bool f;
    double rr;
    for (int i=0;i<vl.vlen;i++) {
        if (vl.v[i]==*l) continue;
        rr=calcos(*l,*r,vl.v[i],dot,&f);
        if (!f) {
            *l=vl.v[i];
            choicedot(vl,vr,dot,l,r);
            return ;
        }
        if (rr==0) {
           if (dcompare(dot[vl.v[i]],dot[*l])) {
               *l=vl.v[i];
               choicedot(vl,vr,dot,l,r);
               return;
           }
        }
    }
    for (int i=0;i<vr.vlen;i++) {
        if (vr.v[i]==*r) continue;
        rr=calcos(*l,*r,vr.v[i],dot,&f);
        if (!f) {
            *r=vr.v[i];
            choicedot(vl,vr,dot,l,r);
            return ;
        }
        if (rr==0) {//角度相同，暂时这样处理，待验证
           if (!dcompare(dot[vr.v[i]],dot[*r])) {
               *r=vr.v[i];
               choicedot(vl,vr,dot,l,r);
               return;
           }
        }
    }
    return;
}

void choicecan(Scandot *rdot,int a1,int a2,SVdot vl,double **dot,int *dl,bool l) {
    int j=0;*dl=-1;
    for (int i=0;i<vl.vlen;i++) {
        if ((vl.v[i]==a1)||(vl.v[i]==a2)) continue;
        bool f;
        double r=calcos(a2,a1,vl.v[i],dot,&f);
        if (f!=l) continue;
        rdot[i].r=r;
        if (*dl==-1) {
            *dl=i;
            rdot[i].less=-1;
            rdot[i].more=-1;
            continue;
        }
        j=*dl;
        while (true) {
            if (r==rdot[j].r) {
                qDebug()<<"samecircle:"<<vl.v[i]<<","<<vl.v[j]<<","<<a1<<","<<a2;
                break;
            }
            if ((r>rdot[j].r)) {
                if (rdot[j].more==-1) {
                    rdot[i].less=j;rdot[i].more=-1;
                    rdot[j].more=i;
                    break;
                }
                j=rdot[j].more;
            } else {
                if (rdot[j].less==-1) {
                    rdot[i].more=j;rdot[i].less=-1;
                    rdot[j].less=i;
                    *dl=i;
                } else {
                    rdot[i].less=j;rdot[i].more=rdot[j].more;
                    rdot[rdot[j].more].less=i;rdot[j].more=i;
                }
                break;
            }
        }
    }    
}

void insertline(SVdot *v,int a1,int a2) {
    if (v->linelen>=v->maxlinelen) {
        v->maxlinelen=v->maxlinelen*2;
        int** temp=new int*[v->maxlinelen];
        memcpy(temp,v->line,v->linelen*sizeof(v->line));
        delete[] v->line;
        v->line=temp;
    }
    v->line[v->linelen]=new int[2];
    if (a1>a2) {
        v->line[v->linelen][0]=a1;
        v->line[v->linelen][1]=a2;
    } else {
        v->line[v->linelen][0]=a2;
        v->line[v->linelen][1]=a1;
    }
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

int isincircle(int a,int b,int c,int d,double** dot) {//d在a,b,c的圆内
    double x0,y0;
    double ax=dot[a][0],bx=dot[b][0],cx=dot[c][0],ay=dot[a][1],by=dot[b][1],cy=dot[c][1];
    x0=( (ax*ax-bx*bx+ay*ay-by*by)*(ay-cy)-(ax*ax-cx*cx+ay*ay-cy*cy)*(ay-by) ) / (2*(ay-cy)*(ax-bx)-2*(ay-by)*(ax-cx));
    y0=( (ax*ax-bx*bx+ay*ay-by*by)*(ax-cx)-(ax*ax-cx*cx+ay*ay-cy*cy)*(ax-bx) ) / (2*(ay-by)*(ax-cx)-2*(ay-cy)*(ax-bx));
    double r=(ax-x0)*(ax-x0) + (ay-y0)*(ay-y0)-((dot[d][0]-x0)*(dot[d][0]-x0)+(dot[d][1]-y0)*(dot[d][1]-y0));
    //r>0,d点在圆心内
    if (r>0) return 1;
    else if (r<0) return -1;
    return 0;
}

SVdot conquer(SVdot *v,SVdot vl,SVdot vr,double** dot) {
    v->linelen=vl.linelen+vr.linelen;
    if (v->linelen>=v->maxlinelen) {
        v->maxlinelen=v->maxlinelen*2;
        int** temp=new int*[v->maxlinelen];
        memcpy(temp,v->line,sizeof (v->line)*(v->linelen));
        delete[] v->line;
        v->line=temp;
    }
    v->linelen=vl.linelen+vr.linelen;
    memcpy(v->line,vl.line,(vl.linelen)*sizeof(vl.line));
    memcpy(v->line+vl.linelen,vr.line,(vr.linelen)*sizeof(vr.line));
    int head,dl=-1,dr=-1;
    int l1=vl.v[vl.vlen-1];
    int r1=vr.v[0];
    choicedot(vl,vr,dot,&l1,&r1);
    Scandot *rdot=new Scandot[vr.vlen],*ldot=new Scandot[vl.vlen];
    while (true) {        
        choicecan(rdot,l1,r1,vr,dot,&head,false);
        int i=head;
        while (i!=-1) {//右边
            if (rdot[i].more==-1) {dr=vr.v[i];break;}
            int r=isincircle(l1,r1,vr.v[i],vr.v[rdot[i].more],dot);
            if (r==-1) {
                dr=vr.v[i];
                break;
            } else {
                delline(v,r1,vr.v[i]);
            }
            i=rdot[i].more;
        }
        if (head==-1) dr=-1;
        choicecan(ldot,r1,l1,vl,dot,&head,true);
        i=head;
        while (i!=-1) {//左边
            if (ldot[i].more==-1) {dl=vl.v[i];break;}
            int r=isincircle(l1,r1,vl.v[i],vl.v[ldot[i].more],dot);
            if (r==-1) {
                dl=vl.v[i];
                break;
            } else {
                delline(v,l1,vl.v[i]);
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
                qDebug()<<"lrcircle:"<<l1<<","<<r1<<","<<vl.v[i]<<","<<vl.v[ldot[i].more];
                break;
            }
        }
        if (dl==-1) r1=dr;
        else l1=dl;
        //下条lr-edge当前是一条直线，就会死循环
    }
    delete[] ldot;
    delete[] rdot;
    delete[] vl.line;
    delete[] vr.line;
    return *v;
}

SVdot divide(SVdot *vdot,double **dot) {
    if (vdot->vlen>3) {
        SVdot vl,vr;
        vl.vlen=vdot->vlen/2;
        vl.v=vdot->v;
        vl.maxlinelen=vl.vlen*2;
        vl.line=new int*[vl.maxlinelen];
        vr.vlen=vdot->vlen-vl.vlen;
        vr.v=vdot->v+vl.vlen;
        vr.maxlinelen=vr.vlen*2;
        vr.line=new int*[vr.maxlinelen];
        return conquer(vdot,divide(&vl,dot),divide(&vr,dot),dot);
    } else {
        if (vdot->vlen==2) {
            vdot->linelen=1;
            vdot->line[0]=new int[2];
            vdot->line[0][0]=vdot->v[0];
            vdot->line[0][1]=vdot->v[1];
        } else {
            vdot->linelen=3;
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
        return SVdot(*vdot);
    }
}

SVdot calDelauney(double **dot,int lens) {
    chaindot chain[lens];
    int dm=0,dl=0,dr=0;
    float dmm=0;

    for (int i=1;i<lens;i++) {
        int j=dm;
        bool f=dcompare(dot[i],dot[j]);
        if (f) {
            while (true) {
                if (chain[j].more==-1) {
                    chain[i].less=j;chain[i].more=-1;
                    chain[j].more=i;
                    dr=i;
                    break;
                }
                if (!dcompare(dot[i],dot[chain[j].more])) {
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
                if (dcompare(dot[i],dot[chain[j].less])) {
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
    vdot.v=new int[lens];
    vdot.vlen=lens;
    vdot.linelen=0;
    vdot.maxlinelen=lens*2;
    vdot.line=new int*[vdot.maxlinelen];
    int i=dl,j=0;
    while (i!=-1) {
        qDebug()<<i<<"x,y:"<<dot[i][0]<<","<<dot[i][1];
        vdot.v[j]=i;
        i=chain[i].more;
        chain[j].sort=i;
        j++;
    }
    qDebug()<<"------------------------";
    return vdot;
}


#endif // DELAUNAY_H
