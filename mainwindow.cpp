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

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    dot=(double**)malloc(len*sizeof(double*));
    for (int i=0;i<len;i++) dot[i]=(double*)malloc(2*sizeof(double));
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_ButtonRam_clicked()
{
    qsrand(QTime::currentTime().msec());
    QPixmap pixmap(ui->labelImg->width(),ui->labelImg->height());
    pixmap.fill(Qt::white);
    QPainter p;
    p.begin(&pixmap);
    p.setPen(Qt::blue);
    QString s;
    for (int i=0;i<len;i++) {
        dot[i][0]=QRandomGenerator::global()->bounded(ui->labelImg->width()-20);
        dot[i][1]=QRandomGenerator::global()->bounded(ui->labelImg->height()-20);
//        dot[i][0]=QRandomGenerator::global()->bounded(len);
//        dot[i][1]=QRandomGenerator::global()->bounded(len);
        p.drawRect(dot[i][0]+10,ui->labelImg->height()-dot[i][1]-10,2,2);
        p.drawText(dot[i][0]+10,ui->labelImg->height()-dot[i][1]-10,QString::number(i));
        s+="("+QString::number(dot[i][0])+","+QString::number(dot[i][1])+"),";
    }
    p.end();
    ui->labelImg->setPixmap(pixmap);
    qDebug()<<s;
}

void MainWindow::on_ButtonCal_clicked()
{
   QElapsedTimer t;
   t.start();
   SVdot vdot=calDelauney(dot,len);
   qDebug()<<"elaspsd:"<<t.elapsed()<<"ms";
//   vdot=divide(&vdot,dot);
   QPixmap pixmap(ui->labelImg->width(),ui->labelImg->height());
   pixmap.fill(Qt::white);
   QPainter p;
   p.begin(&pixmap);
   p.setPen(Qt::red);
   for (int i=0;i<vdot.linelen;i++) {
       if (vdot.line[i]!=nullptr)
       p.drawLine(dot[vdot.line[i][0]][0]+10,ui->labelImg->height()-dot[vdot.line[i][0]][1]-10,dot[vdot.line[i][1]][0]+10,ui->labelImg->height()-dot[vdot.line[i][1]][1]-10);
   }
   p.setPen(Qt::blue);
   for (int i=0;i<len;i++) {
       p.drawRect(dot[i][0]+10,ui->labelImg->height()-dot[i][1]-10,2,2);
       p.drawText(dot[i][0]+10,ui->labelImg->height()-dot[i][1]-10,QString::number(i));
   }
   p.end();
   ui->labelImg->setPixmap(pixmap);
   delvdot(vdot);
}
