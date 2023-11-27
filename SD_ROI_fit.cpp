#include <TCanvas.h>
#include <TGraph.h>
#include <TF1.h>

Double_t SignalByNoiseFunc(Double_t *x, Double_t *par) {
    // 信噪比函数：高斯点源/背景
    // 高斯点源事例数：高斯积分，背景正比于r^2
    // x[0] 是 x 坐标，par[0] \sigma 是高斯点源扩展度，par[1] k 系数
    return x[0]/(par[1] *(1-exp(-x[0]*x[0]/(2*par[0]*par[0]))));
}

void DrawGraphWithFit() {
    // 假设你有两列数据 xData 和 yData
    double xData[] = {1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8};
    double yData[] = {0.4055,0.4074,0.3939,0.3879,0.3879,0.3772,0.3683,0.3682,0.3681,0.3811,0.3759,0.3877,0.3853,0.3933,0.4063,0.4168};
    int numPoints = sizeof(xData) / sizeof(double);

    // 创建一个 TGraph 对象
    TGraph *graph = new TGraph(numPoints, xData, yData);
    graph->SetTitle("Your Graph Title;X Axis Title;Y Axis Title");
    
    // 创建一个自定义的函数
    TF1 *fitFunc = new TF1("fitFunc", SignalByNoiseFunc, 1., 3., 2);
    // 替换 "your_function_formula" 为你自定义的函数表达式，xmin 和 xmax 是拟合范围的最小和最大值
    // 参数初始值
    fitFunc->SetParameter(0, 1.0);
    // fix parameter[0]
    //fitFunc->FixParameter(0, 1.0);
    fitFunc->SetParameter(1, 1.0);

    // 进行数据拟合
    graph->Fit(fitFunc);

    // 创建一个画布并绘制图形
    TCanvas *canvas = new TCanvas("canvas", "Graph with Fit", 800, 600);
    graph->Draw("AP*"); // "AP" 表示绘制图形的同时显示数据点
    fitFunc->Draw("same"); // 在同一画布上绘制拟合函数，使用 "same" 选项保持图形不被覆盖

    // 可以在画布上添加图例等其他绘图选项

    // 保存图形为图像文件（可选）
    //canvas->SaveAs("output.png");
}
