void example() {
    TGraph *gr = new TGraph();
    for (int i = 0; i < 10; i++) {
        gr->SetPoint(i, i, i);
    }
    gStyle->SetOptTitle(1);
    gr->GetXaxis()->SetTitle("X");
    gr->GetYaxis()->SetTitle("Y");
    //gr->GetZaxis()->SetTitle("Z");
    gr->GetXaxis()->SetTitleSize(0.04);
    gr->GetXaxis()->SetTitleFont(42);
    gr->GetXaxis()->SetRangeUser(0,10);
    gr->Draw();
}
