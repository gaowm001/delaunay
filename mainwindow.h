#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QString>
//#include <tetgen.h>
#include <QFileDialog>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

const int len=10000;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:

    void on_ButtonRam_clicked();

    void on_ButtonCal_clicked();

    void on_ButtonTest_clicked();

private:
    Ui::MainWindow *ui;
    QFileInfo filename;
    double** dot;
};
#endif // MAINWINDOW_H
