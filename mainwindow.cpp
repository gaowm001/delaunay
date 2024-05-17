#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QMessageBox>
#include <QApplication>
#include <QRandomGenerator>
#include <QPainter>
#include <QPixmap>
#include <Qtime>
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


void MainWindow::on_ButtonOpen_clicked()
{
    QString file=QFileDialog::getOpenFileName(nullptr,QObject::tr("打开文件"),QDir::current().currentPath(),QObject::tr("poly(*.pole);All File(*)"));
    if (!file.isEmpty()) {
        filename=QFileInfo(file);
        ui->lineEdit->setText(filename.absoluteFilePath());
    }
}

void MainWindow::on_ButtonOut_clicked()
{
    QString sw="pvA";
    filename=ui->lineEdit->text();
    if (ui->checkBoxq->isChecked()) sw+="q";
    if (ui->checkBoxV->isChecked()) sw+=QString::number(ui->SpinBoxV->value())+"V";
    if (ui->checkBoxa->isChecked()) sw+="a"+QString::number(ui->SpinBoxa->value());
    QByteArray bsw=sw.toLocal8Bit(),bfile=filename.absoluteFilePath().toLocal8Bit();
    char* csw=bsw.data(),*file=bfile.data();
    tetgenio in,out;
//    tetgenbehavior *tsw;
//    tsw->parse_commandline(csw);
    if (!in.load_poly(file)) {
        QMessageBox::information(nullptr,"提示","打开文件失败");
        return;
    }
    tetrahedralize(csw,&in,&out);
    out.save_poly("answer");
    out.save_nodes("answer");
    out.save_elements("answer");
    out.save_faces("answer");
    out.save_faces2smesh("answer");
    QMessageBox::information(nullptr,"提示","OK");

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
        dot[i][0]=QRandomGenerator::global()->bounded(120);
        dot[i][1]=QRandomGenerator::global()->bounded(80);
        p.drawRect(dot[i][0]*2,ui->labelImg->height()-dot[i][1]*2,2,2);
        p.drawText(dot[i][0]*2,ui->labelImg->height()-dot[i][1]*2,QString::number(i));
        s+="("+QString::number(dot[i][0])+","+QString::number(dot[i][1])+"),";
    }
    p.end();
    ui->labelImg->setPixmap(pixmap);
    qDebug()<<s;
}

void MainWindow::on_ButtonCal_clicked()
{
   SVdot vdot=calDelauney(dot,len);
   vdot=divide(&vdot,dot);
   QPixmap pixmap(ui->labelImg->width(),ui->labelImg->height());
   pixmap.fill(Qt::white);
   QPainter p;
   p.begin(&pixmap);
   p.setPen(Qt::blue);
   for (int i=0;i<len;i++) {
       p.drawRect(dot[i][0]*2,ui->labelImg->height()-dot[i][1]*2,2,2);
       p.drawText(dot[i][0]*2,ui->labelImg->height()-dot[i][1]*2,QString::number(i));
   }
   p.setPen(Qt::red);
   for (int i=0;i<vdot.linelen;i++) {
       if (vdot.line[i]!=nullptr)
       p.drawLine(dot[vdot.line[i][0]][0]*2,ui->labelImg->height()-dot[vdot.line[i][0]][1]*2,dot[vdot.line[i][1]][0]*2,ui->labelImg->height()-dot[vdot.line[i][1]][1]*2);
   }
   p.end();
   ui->labelImg->setPixmap(pixmap);

}
