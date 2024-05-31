#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QMessageBox>
#include <QApplication>
#include <QRandomGenerator>
#include <QPainter>
#include <QPixmap>
#include <Qtime>
#include <QElapsedTimer>
#include "delaunay.h"
#include <QDebug>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    dot=(double**)malloc(len*sizeof(double*));
    for (int i=0;i<len;i++) dot[i]=(double*)malloc(2*sizeof(double));
    QPixmap pixmap(ui->labelImg->width(),ui->labelImg->height());
    pixmap.fill(Qt::white);
    ui->labelImg->setPixmap(pixmap);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_ButtonRam_clicked()
{
//    dot[0][0]=-1;dot[0][1]=-1;
//    dot[1][0]=-1;dot[1][1]=1;
//    dot[2][0]=0;dot[2][1]=0;
//    dot[3][0]=1;dot[3][1]=-1;
//    dot[4][0]=1;dot[4][1]=1;
//    SVdot vdot=calDelauney(dot,len);
//    for (int i=0;i<vdot.linelen;i++) delete[] vdot.line[i];
//    delete[] vdot.line;
//    vdot.linelen=4;
//    vdot.maxlinelen=vdot.linelen*2;
//    vdot.line=new int*[vdot.maxlinelen];
//    vdot.line[0]=new int[2];vdot.line[0][0]=0;vdot.line[0][1]=1;
//    vdot.line[1]=new int[2];vdot.line[0][0]=1;vdot.line[0][1]=3;
//    vdot.line[2]=new int[2];vdot.line[0][0]=3;vdot.line[0][1]=4;
//    vdot.line[3]=new int[2];vdot.line[0][0]=4;vdot.line[0][1]=5;


    qsrand(QTime::currentTime().msec());
    QPixmap pixmap(ui->labelImg->width(),ui->labelImg->height());
    pixmap.fill(Qt::white);
    QPainter p;
    p.begin(&pixmap);
    p.setPen(Qt::blue);
    QString s;
    for (int i=0;i<len;i++) {
        if (len<ui->labelImg->width()-20) {
            dot[i][0]=QRandomGenerator::global()->bounded(ui->labelImg->width()-20);
            dot[i][1]=QRandomGenerator::global()->bounded(ui->labelImg->height()-20);
        } else {
            dot[i][0]=QRandomGenerator::global()->bounded(len);
            dot[i][1]=QRandomGenerator::global()->bounded(len);
        }
        p.drawRect(dot[i][0]+10,ui->labelImg->height()-dot[i][1]-10,2,2);
        p.drawText(dot[i][0]+10,ui->labelImg->height()-dot[i][1]-10,QString::number(i));
        s+="("+QString::number(dot[i][0])+","+QString::number(dot[i][1])+"),";
    }
    p.end();
    ui->labelImg->setPixmap(pixmap);
    ui->ButtonCal->setEnabled(true);
    qDebug()<<s;
}

void MainWindow::on_ButtonCal_clicked()
{
   QElapsedTimer t;
   t.start();
   SVdot vdot=calDelaunay(dot);
   ui->textEdit->append("elasped:"+QString::number(t.elapsed())+"ms");
   t.invalidate();
   QPixmap pixmap(ui->labelImg->width(),ui->labelImg->height());
   pixmap.fill(Qt::white);
   QPainter p;
   p.begin(&pixmap);
   p.setPen(Qt::red);
   for (int i=0;i<vdot.linelen;i++) {
       p.drawLine(dot[vdot.line[i][0]][0]+10,ui->labelImg->height()-dot[vdot.line[i][0]][1]-10,dot[vdot.line[i][1]][0]+10,ui->labelImg->height()-dot[vdot.line[i][1]][1]-10);
   }
   p.setPen(Qt::blue);
   for (int i=0;i<len;i++) {
       p.drawRect(dot[i][0]+10,ui->labelImg->height()-dot[i][1]-10,2,2);
       p.drawText(dot[i][0]+10,ui->labelImg->height()-dot[i][1]-10,QString::number(i));
   }
   p.end();
   ui->labelImg->setPixmap(pixmap);
   qApp->processEvents();
   if (checkDelaunay1(vdot,dot)) ui->textEdit->append("All hollow circle") ;//判断空圆
   else QMessageBox::information(nullptr,"Info","Delauney is checked error!");
   qApp->processEvents();
   if (checkDelaunay(vdot,dot)) ui->textEdit->append("No intersect");//判断相交
   else QMessageBox::information(nullptr,"Info","Delauney is checked error!");
   delvdot(vdot);
}

void MainWindow::on_ButtonTest_clicked()
{
    int count=100,j1=0,j2=0,k=0;
    int *test=new int[count];
    QString s1;
    for (int i=0;i<count;i++) {
        QString s;
    //        int adot[10][2]={17,77,242,266,221,18,115,177,0,94,130,174,199,191,100,24,260,273,0,223};
    //        for (int j=0;j<len;j++) {dot[j][0]=adot[j][0];dot[j][1]=adot[j][1];}
        for (int j=0;j<len;j++) {
            dot[j][0]=QRandomGenerator::global()->bounded(len);
            dot[j][1]=QRandomGenerator::global()->bounded(len);
            s+="("+QString::number(dot[j][0])+","+QString::number(dot[j][1])+"),";
        }
        QElapsedTimer t;
        t.start();
        SVdot vdot=calDelaunay(dot);
        k+=t.elapsed();
        t.invalidate();
        test[i]=0;
/*
        if (!checkDelaunay1(vdot,dot)) {
            QPixmap pixmap(ui->labelImg->width(),ui->labelImg->height());
            pixmap.fill(Qt::white);
            QPainter p;
            p.begin(&pixmap);
            p.setPen(Qt::red);
            for (int i=0;i<vdot.linelen;i++) {
                p.drawLine(dot[vdot.line[i][0]][0]*(ui->labelImg->width()-20)/len+10,ui->labelImg->height()-dot[vdot.line[i][0]][1]*(ui->labelImg->height()-20)/len-10,dot[vdot.line[i][1]][0]*(ui->labelImg->width()-20)/len+10,ui->labelImg->height()-dot[vdot.line[i][1]][1]*(ui->labelImg->height()-20)/len-10);
            }
            p.setPen(Qt::blue);
            for (int i=0;i<len;i++) {
                p.drawRect(dot[i][0]*(ui->labelImg->width()-20)/len+10,ui->labelImg->height()-dot[i][1]*(ui->labelImg->height()-20)/len-10,2,2);
                p.drawText(dot[i][0]*(ui->labelImg->width()-20)/len+10,ui->labelImg->height()-dot[i][1]*(ui->labelImg->height()-20)/len-10,QString::number(i));
            }
            p.end();
            ui->labelImg->setPixmap(pixmap);
            qDebug()<<s;
            s1=s;
            qApp->processEvents();
            delvdot(vdot);

            vdot=calDelaunay(dot);
            test[i]=1;
            j1++;
            qDebug()<<QString::number(i)<<":circle error";
        }
        if (!checkDelaunay(vdot,dot)) {test[i]=2;j2++;ui->textEdit->append(QString::number(i)+":line error");
        }
*/
        delvdot(vdot);
        ui->textEdit->append(QString::number(i)+" end");
        qApp->processEvents();
    }    
    qDebug()<<s1;
    ui->textEdit->append("len:"+QString::number(len)+" total:"+QString::number(count)+" circle error:"+QString::number(j1)+" line error:"+QString::number(j2)+" average time(ms):"+QString::number(k/count));
    delete[] test;
}
