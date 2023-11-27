double TDCcorrect(double pe){
    double dt,drt;
    //if(pe>200)return -5.880;
    if(pe>500)return -7.37;
    dt=20*pow(pe,-0.1)+10*pow(pe,-0.08);
    //drt=20*pow(10.,-0.1)+10*pow(10,-0.08);   //24.204202
    return dt-24.20;
}
