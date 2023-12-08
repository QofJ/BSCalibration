#include "BSTimeoffCalibration.h"
#include "astro-csz.h"
#define PETH 0.2

BSTimeoffCalibration *BSTimeoffCalibration::m_myself = 0;
G4KM2A_Geometry *BSTimeoffCalibration::geom = 0;
LHEvent *BSTimeoffCalibration::lhevent = 0;
int BSTimeoffCalibration::ITER = 0;
TH2D *BSTimeoffCalibration::h2 = new TH2D("h2", "time residual of ED", 1000, -50, 50, 5500, 0, 5500);

BSTimeoffCalibration::BSTimeoffCalibration()
{
    geom = G4KM2A_Geometry::GetInstance(7);
    Init();
}

int BSTimeoffCalibration::Init()
{
    int i;
    double ra, dec;
    for (i = 0; i < MNED; i++)
    {
        conical_residualsED[i] = 0;
        if (ITER == 0)
            ToffsetED[i] = 0;
        cal_statED[i] = 0;
    }
    theta = -1;
    phi = -1;
    rchi = 0;
    rchi_sum = 0;
    dir[0] = 0;
    dir[1] = 0;
    dir[2] = 1;
    dir_set[0] = -10;
    dir_set[1] = -10;
    dir_set[2] = -10;
    core[0] = 0;
    core[1] = 0;
    a = 0;
    c0 = 0;
    Sigma = 0;
    Sigma_mean = 0;
    Nevents = 0;
    Dazcos[0] = 0;
    Dazcos[1] = 0;
    Dazcos[2] = 0;

    BSTimeoffCalibration::h2->Reset();

    return 0;
}

int BSTimeoffCalibration::SetEvent(LHEvent *lhe)
{
    lhevent = lhe;
    return 0;
}

int BSTimeoffCalibration::SetITER(int iter)
{
    ITER = iter;
    return 0;
}

int BSTimeoffCalibration::SetToffset(const char *file)
{
    int id;
    double toff;
    char type[20];
    FILE *f = fopen(file, "r");
    if (f == NULL)
    {
        std::cout << "Error opening file!" << std::endl;
        return -1;
    }
    else
    {
        std::cout << "File opened successfully!" << std::endl;
    }
    fscanf(f, "%*[^\n]%*c");

    for (int i = 0; i < MNED; i++)
    {
        ToffsetED[i] = 0;
    }
    while (!feof(f))
    {
        fscanf(f, "%s %d %lf\n", type, &id, &toff);
        if (strcmp(type, "ED") != 0)
            break;
        if (id < 0 || id >= MNED)
            continue;
        ToffsetED[id] = toff;
    }
    fclose(f);
    return 0;
}

int BSTimeoffCalibration::SetC0(double C0)
{
    this->c0 = C0;
    return 0;
}

int BSTimeoffCalibration::SetSource(double mjd, double ra, double dec, int flag)
{
    double zen, az;
    double Mjd;
    if (mjd < 0)
        Mjd = lhevent->GetMjd();
    else
        Mjd = mjd;
    precession2New(Mjd, &ra, &dec);
    e2h_argo(Mjd, ra, dec, &zen, &az);
    theta = zen * DD2R;
    phi = az * DD2R;
    if (flag == 0)
    {
        dir_set[0] = sin(theta) * cos(phi);
        dir_set[1] = sin(theta) * sin(phi);
        dir_set[2] = cos(theta);
        this->c0_set = this->c0;
    }
    else if (flag == 1)
    {
        dir_set[0] = sin(theta) * cos(phi);
        dir_set[1] = sin(theta) * sin(phi);
        dir_set[2] = cos(theta);
        this->c0_set = this->c0 - (dir_set[0] - dir[0]) * (this->core[0]) - (dir_set[1] - dir[1]) * (this->core[1]);
    }
    else
    {
        dir[0] = sin(theta) * cos(phi);
        dir[1] = sin(theta) * sin(phi);
        dir[2] = cos(theta);
    }
    return 0;
}

int BSTimeoffCalibration::SetCore(double core_x, double core_y)
{
    core[0] = core_x;
    core[1] = core_y;
    return 0;
}

int BSTimeoffCalibration::Setnparticle(TClonesArray &tHits, int np, double pe)
{
    LHHit *tHit;
    int i;
    double ne;
    for (i = 0; i < np; i++)
    {
        tHit = (LHHit *)((tHits)[i]);
        if (tHit->GetStatus() > -1)
        {
            ne = tHit->GetPe() / pe;
            if (ne < PETH)
                ne = 0;
            tHit->SetPe(ne);
        }
    }
    return 0;
}

int BSTimeoffCalibration::SetBadHit(TClonesArray &tHits, int np, const char *style)
{
    int i, id;
    LHHit *tHit;
    if (strcmp(style, "ED") == 0)
    {
        for (i = 0; i < np; i++)
        {
            tHit = (LHHit *)((tHits)[i]);
            id = tHit->GetId();
            if (geom->GetEDId2(id) < 0)
                tHit->SetStatus(-2);
        }
    }
    else if (strcmp(style, "MD") == 0)
    {
        for (i = 0; i < np; i++)
        {
            tHit = (LHHit *)((tHits)[i]);
            id = tHit->GetId();
            if (geom->GetMDId2(id) < 0)
                tHit->SetStatus(-2);
        }
    }
    return 0;
}

double BSTimeoffCalibration::GetTheta()
{
    return this->theta;
}

int BSTimeoffCalibration::SetPreToffset(const char *meanCalfile, const char *Pre_Toffsetfile)
{
    int id;
    double time;
    FILE *f1 = fopen(meanCalfile, "r");
    FILE *f2 = fopen(Pre_Toffsetfile, "r");

    fscanf(f1, "%*[^\n]%*c");
    fscanf(f2, "%*[^\n]%*c");

    for (int i = 0; i < MNED; i++)
    {
        Pre_ToffsetED[i] = -1000;
    }

    while (!feof(f2))
    {
        fscanf(f2, "ED %d %lf\n", &id, &time);
        Pre_ToffsetED[id] = 0 - time;
    }
    fclose(f2);

    while (!feof(f1))
    {
        fscanf(f1, "ED %d %lf\n", &id, &time);
        if (Pre_ToffsetED[id] > -500)
            Pre_ToffsetED[id] += time;
    }

    return 0;
}

double BSTimeoffCalibration::GetPhi()
{
    return this->phi;
}

double BSTimeoffCalibration::GetA()
{
    return this->a;
}

double BSTimeoffCalibration::GetC0()
{
    return this->c0;
}

double BSTimeoffCalibration::GetSigma()
{
    return this->Sigma;
}

double BSTimeoffCalibration::GetC0_set()
{
    return this->c0_set;
}

double *BSTimeoffCalibration::GetstatED()
{
    return cal_statED;
}

double *BSTimeoffCalibration::GetToffsetED()
{
    return ToffsetED;
}

int BSTimeoffCalibration::Addresidual(TClonesArray &tHits, int np, const char *style)
{
    int i, nus = 0;
    double dirl, dirm, dirn, dr, dt, tt = 0;
    double drc, dtc, w_square;
    double x, y, z, dx, dy, dz;
    LHHit *tHit;
    dirl = dir_set[0];
    dirm = dir_set[1];
    dirn = dir_set[2];

    /*double dtt[np];
    int dtt_ID[np] = { 0 };
    double mean = 0;
    double w_con = 0;*/

    // double reftime = GetRefTime(1);

    for (i = 0; i < np; i++)
    {
        tHit = (LHHit *)((tHits)[i]);
        if (tHit->GetStatus() < 2.5 || tHit->GetPe() < 1.e-3)
            continue;
        // printf("\n%d %d ",i,tHit->GetId());
        if (geom->Getxyz(tHit->GetId(), x, y, z, 1, style) < 0)
            continue;
        w_square = 1;
        dx = x - core[0];
        dy = y - core[1];
        dz = z;
        dr = sqrt(dx * dx + dy * dy + dz * dz - pow(dx * dirl + dy * dirm - dz * dirn, 2.));

        // if (dr > 100)continue;

        dt = tHit->GetTime() - ToffsetED[tHit->GetId()] - tt - (dirl * x + dirm * y - dirn * z + a * dr + c0_set) / 0.2998;

        // printf("ED %d %lf\n", tHit->GetId(), dt);

        // dt -= reftime;

        // dt = ((dir[0] - dirl) * dx + (dir[1] - dirm) * dy - (dir[2] - dirn) * dz) / 0.2998;

        // dt = tHit->GetTime() - tt - (dirl * x + dirm * y - dirn * z + a * dr + c0_set) / 0.2998
        //     + ((dir[0] - dirl) * dx + (dir[1] - dirm) * dy - (dir[2] - dirn) * dz) / 0.2998;

        /*if (ITER >= 0)
        {
            if (fabs(dt - ToffsetED[tHit->GetId()]) > 10)continue;
        }*/

        // if (dt > 0 || dt < -20) w_square = w_square * exp(-1.0 * pow(dt / 10., 2));

        /*if (ITER > 0)
        {
            drc = sqrt(dx * dx + dy * dy + dz * dz - pow(dx * dir[0] + dy * dir[1] - dz * dir[2], 2.));
            //dtc = dt - ToffsetED[tHit->GetId()];
            dtc = tHit->GetTime() - ToffsetED[tHit->GetId()] - tt - (dir[0] * x + dir[1] * y - dir[2] * z + a * drc + c0) / 0.2998;
            if (dtc > 0 || dtc < -20) w_square = w_square * exp(-1.0 * pow(dtc / 10., 2));

            //if (dtc > 5 || dtc < -5)continue;
        }*/

        /* dtt[nus] = dt;
         dtt_ID[nus] = tHit->GetId();
         mean += w_square * ToffsetED[tHit->GetId()];
         w_con += w_square;*/
        nus++;
        // dt = tHit->GetTime()-tt-(dirl*x+dirm*y-dirn*z+trec->rec_c0)/0.2998;
        // w_square = w_square / (this->Sigma * (this->NfiltED - 1));
        cal_statED[tHit->GetId()] += w_square;
        // if (delta != 0)dt -= ToffsetED[tHit->GetId()];
        // conical_residualsED[tHit->GetId()] += w_square * dt;

        BSTimeoffCalibration::h2->Fill(dt, tHit->GetId() - 0.5);

        // if (dt < -15)printf("dt = %lf, sigma = %lf\n", dt, this->Sigma);
    }

    // printf("Over\n");

    /*if (nus == 0)return -1;
    mean /= w_con;
    for (i = 0; i < nus; i++)
    {
        BSTimeoffCalibration::h2->Fill(dtt[i] + mean, dtt_ID[i] - 0.5);
    }*/

    Nevents++;
    Dazcos[0] += dir[0] - dir_set[0];
    Dazcos[1] += dir[1] - dir_set[1];
    Dazcos[2] += dir[2] - dir_set[2];

    this->Sigma_mean += (this->Sigma - this->Sigma_mean) / Nevents;

    return nus;
}

int BSTimeoffCalibration::GenerateToffset(const char *style, int MPmode, double delta, double step)
{
    printf("Sigma_Mean = %lf\n", this->Sigma_mean);

    if (strcmp(style, "ED") == 0)
    {
        if (MPmode == 0)
        {

            for (int i = 0; i < MNED; i++)
            {
                if (cal_statED[i] == 0)
                    ToffsetED[i] = 0;
                else
                    ToffsetED[i] += step * (conical_residualsED[i] / cal_statED[i]);
            }
        }
        else if (MPmode == 1)
        {
            for (int i = 0; i < MNED; i++)
            {
                if (cal_statED[i] == 0)
                    continue;
                // else ToffsetED[i] = GetMPValue(i, delta);
                else
                    ToffsetED[i] += step * (GetMPValue(i, delta) - ToffsetED[i]);
            }
        }
        else
        {
            for (int i = 0; i < MNED; i++)
            {
                if (cal_statED[i] == 0)
                    continue;
                // else ToffsetED[i] = GetMPValue(i, delta);
                else
                    ToffsetED[i] += step * GetMPValue(i, delta);
            }
            // ApplyDaz();
        }
    }
    return 0;
}

double BSTimeoffCalibration::NormToffset()
{
    int i;
    int num = 0;
    double mean = 0;
    // if (ITER != 0)
    //{
    for (i = 0; i < MNED; i++)
    {
        if (cal_statED[i] == 0)
            continue;
        num++;
        mean += ToffsetED[i];
    }
    if (num != 0)
    {
        mean /= num;
        for (i = 0; i < MNED; i++)
        {
            if (cal_statED[i] == 0)
                continue;
            ToffsetED[i] -= mean;
        }
    }
    //}
    return mean;
}

int BSTimeoffCalibration::spacetimefilter(TClonesArray &tHits, int np, int twind, int rwind, const char *style)
{
    int i, j, t, id, tflag, nflag, maxid, maxt;
    double a, b, x, y, z, r;
    LHHit *tHit;
    int Nbin, Nstep, Nstar, Nall = 0;
    double flag, rmax, tdc[10000];
    if (strcmp(style, "WCDA") == 0)
    {
        Nbin = 2200 - 400;
        Nstep = 10;
        Nstar = 600;
    }
    else
    {
        Nbin = 10000 - 1000;
        Nstep = 50;
        Nstar = 4000;
        // if(np> 5*(geom->GetNED())){
        // Nbin=10000-4600;
        //   rwind=1000;
        //}
    }
    for (i = 0; i < Nbin; i++)
        tdc[i] = 0;
    for (i = 0; i < np; i++)
    {
        tHit = (LHHit *)((tHits)[i]);
        if (tHit->GetStatus() > -1)
        {
            if (tHit->GetPe() < 1.e-3)
                tHit->SetStatus(0);
            else
            {
                tHit->SetStatus(5);
                t = int(tHit->GetTime() - ToffsetED[tHit->GetId()]);
                if (t >= 0 && t < Nbin)
                {
                    tdc[t]++;
                    Nall++;
                }
                // if(t>=0&&t<Nbin){ tdc[t]+=tHit->GetPe(); Nall++;}
            }
        }
    }
    if (Nall < 4)
        return 0;
    // to search for the maximum circle with radius rwind
    rmax = 0;
    for (tflag = Nstar; tflag < Nbin - twind; tflag += Nstep)
    {
        while (tdc[tflag] == 0 && tflag < Nbin - twind - 1)
            tflag++;
        // if(rmax>3){
        //    flag=0;
        //    for(j=tflag;j<tflag+twind;j++)if(j<Nbin-twind)flag +=tdc[j];
        //    if(flag<rmax)continue;
        // }
        flag = 0;
        nflag = 0;
        for (i = 0; i < np; i++)
        {
            tHit = (LHHit *)((tHits)[i]);
            if (tHit->GetStatus() > 0 && tHit->GetTime() - ToffsetED[tHit->GetId()] >= tflag && tHit->GetTime() - ToffsetED[tHit->GetId()] < tflag + twind)
            {
                if (geom->Getxyz(tHit->GetId(), x, y, z, 1, style) < 0)
                    continue;
                a = x;
                b = y;
                id = tHit->GetId();
                flag = 0;
                nflag = 0;
                for (j = 0; j < np; j++)
                {
                    tHit = (LHHit *)((tHits)[j]);
                    if (tHit->GetStatus() > 0 && tHit->GetTime() - ToffsetED[tHit->GetId()] >= tflag && tHit->GetTime() - ToffsetED[tHit->GetId()] < tflag + twind)
                    {
                        if (geom->Getxyz(tHit->GetId(), x, y, z, 1, style) < 0)
                            continue;
                        r = sqrt(pow(x - a, 2) + pow(y - b, 2));
                        if (r < rwind)
                        {
                            flag += tHit->GetPe();
                            nflag++;
                        }
                    }
                }
                if (flag > rmax && nflag > 10)
                {
                    rmax = flag;
                    maxid = id;
                    maxt = tflag;
                }
                if (nflag > rmax)
                {
                    rmax = nflag;
                    maxid = id;
                    maxt = tflag;
                }
            }
        } // end of np
    }     // end of tflag
    // filter out noise
    /*
    if(rmax<4){
            printf("ramx %.2lf %d %d rwind=%d twin=%d \n",rmax,maxt,maxid,rwind,twind);
            for(j=0;j<np;j++){
                  tHit = (LHHit *)((tHits)[j]);
                  printf("   %d %.2lf %.4lf %d\n",tHit->GetId(),tHit->GetTime(),tHit->GetPe(),tHit->GetStatus());
                }
     }
    */
    if (maxid > 0 && geom->Getxyz(maxid, x, y, z, 1, style) > 0)
    {
        a = x;
        b = y;
        for (i = 0; i < np; i++)
        {
            tHit = (LHHit *)((tHits)[i]);
            if (tHit->GetStatus() > 0)
            {
                if (tHit->GetTime() - ToffsetED[tHit->GetId()] < maxt || tHit->GetTime() - ToffsetED[tHit->GetId()] >= maxt + twind)
                    tHit->SetStatus(0);
                else
                {
                    if (geom->Getxyz(tHit->GetId(), x, y, z, 1, style) > 0)
                    {
                        r = sqrt(pow(x - a, 2) + pow(y - b, 2));
                        if (r > rwind)
                            tHit->SetStatus(0);
                    }
                }
            }
        }
    }
    // return maxt;
    return (int)rmax;
}

int BSTimeoffCalibration::SetA(double alpha)
{
    this->a = alpha;
    return 0;
}

int BSTimeoffCalibration::noisefilter(TClonesArray &tHits, int np, int twindl, int twindh, int rwind, const char *style)
{
    int i, nus = 0;
    double dt, dirl, dirm, dirn, dr, tt = 0;
    double x, y, z, dx, dy, dz;
    LHHit *tHit;
    dirl = dir[0];
    dirm = dir[1];
    dirn = dir[2];
    for (i = 0; i < np; i++)
    {
        tHit = (LHHit *)((tHits)[i]);
        if (tHit->GetStatus() < 0 || tHit->GetPe() < 1.e-3)
            continue;
        if (geom->Getxyz(tHit->GetId(), x, y, z, 1, style) < 0)
            continue;
        dx = x - core[0];
        dy = y - core[1];
        dz = z;
        dr = sqrt(dx * dx + dy * dy + dz * dz - pow(dx * dirl + dy * dirm - dz * dirn, 2.));
        dt = tHit->GetTime() - ToffsetED[tHit->GetId()] - tt - (dirl * x + dirm * y - dirn * z + a * dr + c0) / 0.2998;

        if (dt > twindl && dt < twindh)
        {
            if (dr < rwind)
            {
                tHit->SetStatus(5);
                nus++;
            }
            else
                tHit->SetStatus(2);
        }
        else
        {
            if (dr < rwind)
                tHit->SetStatus(1);
            else
                tHit->SetStatus(0);
        }
    }
    return nus;
}

int BSTimeoffCalibration::conicalfit(TClonesArray &tHits, int np, float alpha, const char *style, int IsSet)
{
    int i, j, k, nus, AA, Flag = 0;
    int BinaryFlag = 0;
    double dt, w, c[2];
    double x, y, z, dx, dy, dz, dr, sigma;
    double dirl, dirm, dirn;
    double a_temp, c0_temp = 0;
    double c0_u = -1, c0_d = -1;
    if (alpha > 0)
        AA = 1;
    else
        AA = 0;
    TMatrixD A(2, 2);
    TVectorD b(2), r(2);
    LHHit *tHit;

    if (IsSet == 0)
    {
        dirl = dir[0];
        dirm = dir[1];
        dirn = dir[2];
    }
    else
    {
        dirl = dir_set[0];
        dirm = dir_set[1];
        dirn = dir_set[2];
    }

    A.Zero();
    b.Zero();
    nus = 0;
    // for (int iter = 0; iter < 200; iter++)
    for (int iter = 0; iter < 1000; iter++)
    {
        nus = 0;
        sigma = 0;
        for (i = 0; i < np; i++)
        {
            tHit = (LHHit *)((tHits)[i]);
            if (tHit->GetStatus() < 2.5 || tHit->GetPe() < 1.e-3)
                continue;
            if (geom->Getxyz(tHit->GetId(), x, y, z, 1, style) < 0)
                continue;
            dx = x - core[0];
            dy = y - core[1];
            dz = z;
            dr = (dx * dx + dy * dy + dz * dz - pow(dx * dirl + dy * dirm - dz * dirn, 2.));
            if (dr < 0)
                dr = 0;
            dr = sqrt(dr);

            w = 1.;
            // if (ITER > 0 && iter > 0)
            if (iter > 0)
            {
                // dt = tHit->GetTime() - ToffsetED[tHit->GetId()] - (dir[0] * x + dir[1] * y - dir[2] * z + c0) / 0.2998;
                dt = tHit->GetTime() - ToffsetED[tHit->GetId()] - (dirl * x + dirm * y - dirn * z + a_temp * dr + c0_temp) / 0.2998;

                if (dt > 0 || dt < -20)
                    w = w * exp(-0.5 * pow(dt / 10., 2));

                sigma += dt * dt;
            }
            nus++;

            c[0] = dr * (1 - AA) * w;
            c[1] = 1 * w;
            for (k = 0; k < 2; k++)
            {
                for (j = 0; j < 2; j++)
                    A[j][k] += c[j] * c[k];
                b[k] += ((tHit->GetTime() - ToffsetED[tHit->GetId()]) * 0.2998 - AA * alpha * dr - dirl * x - dirm * y + dirn * z) * w * c[k];
            }
        }

        // if (iter == 1)
        //{
        //     sigma = sqrt(sigma / (nus - 1.));
        //     break;
        // }

        if (nus < 2)
            break;
        TDecompSVD svd(A);
        Bool_t ok;
        r = svd.Solve(b, ok);

        if (!ok)
            break;

        if (ok)
        {
            // printf("C0 = %lf\n", c0_temp);

            this->NfiltED = nus;
            Flag = 1;
            sigma = sqrt(sigma / (nus - 1.));
            if (AA == 0)
            {
                a_temp = r[0];
                if (fabs(r[1] - c0_temp) < 0.0001)
                {
                    c0_temp = r[1];
                    break;
                }
                else
                {
                    c0_temp = r[1];
                }
            }
            else
            {
                a_temp = alpha;

                if (fabs(r[1] - c0_temp) < 0.0001)
                {
                    c0_temp = r[1];
                    break;
                }
                else
                {
                    c0_temp = r[1];
                }

                /*if (fabs(r[1] - c0_temp) > 0.1 && BinaryFlag == 0)
                {
                    c0_temp = r[1];
                }
                else if (BinaryFlag == 0)
                {
                    BinaryFlag = 1;
                }
                if (BinaryFlag == 1)
                {
                    if (c0_u < 0 && c0_d < 0)
                    {
                        if (r[1] - c0_temp < 0)
                        {
                            c0_u = c0_temp;
                            c0_temp -= 1;
                        }
                        else
                        {
                            c0_d = c0_temp;
                            c0_temp += 1;
                        }
                    }
                    else if (c0_u > 0 && c0_d < 0)
                    {
                        if (r[1] - c0_temp < 0)
                        {
                            c0_u = c0_temp;
                            c0_temp -= 1;
                        }
                        else
                        {
                            c0_d = c0_temp;
                            c0_temp = (c0_u + c0_d) / 2;
                            BinaryFlag = 2;
                        }
                    }
                    else if (c0_u < 0 && c0_d>0)
                    {
                        if (r[1] - c0_temp > 0)
                        {
                            c0_d = c0_temp;
                            c0_temp += 1;
                        }
                        else
                        {
                            c0_u = c0_temp;
                            c0_temp = (c0_u + c0_d) / 2;
                            BinaryFlag = 2;
                        }
                    }
                }
                if (BinaryFlag == 2)
                {
                    if (r[1] - c0_temp < 0)
                    {
                        c0_u = c0_temp;
                        c0_temp = (c0_u + c0_d) / 2;
                    }
                    else
                    {
                        c0_d = c0_temp;
                        c0_temp = (c0_u + c0_d) / 2;
                    }
                    if (fabs(c0_u - c0_d) < 1e-9)break;
                }*/
            }
        }
        // if (ITER == 0)break;
    }
    // printf("Over\n");

    if (Flag == 1)
    {
        this->Sigma = sigma;

        if (IsSet == 0)
        {
            this->a = a_temp;
            this->c0 = c0_temp;
        }
        else
        {
            this->a = a_temp;
            this->c0_set = c0_temp;
        }
        return 0;
    }
    else
        return -1;
}

TH2D *BSTimeoffCalibration::GetHist()
{
    return BSTimeoffCalibration::h2;
}

int BSTimeoffCalibration::conicalfit_dir(TClonesArray &tHits, int np, float alpha, const char *style, int nmode, double mjd, double ra0, double dec0)
{
    double ra, dec;

    if (np < 4)
        return -1;
    int i, j, k, nus, AA, Flag = 0;
    double dt, w, the, phi_temp, dirl, dirm, dirn, dr, c[4];
    double x, y, z, dx, dy, dz, sigma;
    double a_temp = alpha;
    if (alpha > 0)
        AA = 1;
    else
        AA = 0; // to fit the alpha
    TMatrixD A(4, 4);
    TVectorD b(4), r(4);
    LHHit *tHit;

    if (nmode != 0)
    {
        precession2New(mjd, &ra0, &dec0);
        e2h_argo(mjd, ra0, dec0, &the, &phi_temp);
        the = the * DD2R;
        phi_temp = phi_temp * DD2R;
        dirl = sin(the) * cos(phi_temp);
        dirm = sin(the) * sin(phi_temp);
        dirn = cos(the);
    }
    else
    {
        dirl = dir[0];
        dirm = dir[1];
        dirn = dir[2];
    }

    for (int iter = 0; iter < 30; iter++)
    {
        for (i = 0; i < 4; i++)
        {
            for (j = 0; j < 4; j++)
                A[i][j] = 0;
            b[i] = 0;
        }
        nus = 0;
        sigma = 0;
        for (i = 0; i < np; i++)
        {
            tHit = (LHHit *)((tHits)[i]);
            if (tHit->GetStatus() < 2.5 || tHit->GetPe() < 1.e-3)
                continue; // not use the noise hit filter by planar fit
            if (geom->Getxyz(tHit->GetId(), x, y, z, 1, style) < 0)
                continue;
            dx = x - core[0];
            dy = y - core[1];
            dz = z;
            dr = (dx * dx + dy * dy + dz * dz - pow(dx * dirl + dy * dirm - dz * dirn, 2.));
            if (dr < 0)
                dr = 0;
            dr = sqrt(dr);
            w = 1.;
            if (iter > 0)
            {
                // dt = tHit->GetTime()-(r[0]*x+r[1]*y-dirn*z+(r[2]+alpha*AA)*dr+r[3])/0.2998; //ns
                // dt = tHit->GetTime() - ToffsetED[tHit->GetId()] - (r[0] * x + r[1] * y - dirn * z + r[3]) / 0.2998; //ns please use this dt
                dt = tHit->GetTime() - ToffsetED[tHit->GetId()] - (r[0] * x + r[1] * y - dirn * z + a_temp * dr + r[3]) / 0.2998;
                if (dt > 0 || dt < -20)
                    w = w * exp(-0.5 * pow(dt / 10., 2)); // sig=10ns is better than 5,15,20
                // printf("id = %d, dt = %lf, w = %lf\n", tHit->GetId(), dt, w);
                sigma += dt * dt;
            }
            nus++;
            // w=sqrt(w); has a little worse for this
            // w*=sqrt(tHit->GetNp()); not good
            c[0] = x * w;
            c[1] = y * w;
            c[2] = dr * (1 - AA) * w;
            c[3] = 1 * w;
            for (k = 0; k < 4; k++)
            {
                for (j = 0; j < 4; j++)
                    A[j][k] += c[j] * c[k];
                b[k] += ((tHit->GetTime() - ToffsetED[tHit->GetId()]) * 0.2998 - AA * alpha * dr + dirn * z) * w * c[k];
            }
            // if(DEBUG)printf("%d %lf %lf %lf %lf,%lf\n",i,c[0],c[1],tHit->GetTime(),dr,dirn*z);
        }
        // solve the equation
        if (nus < 4)
            break;
        TDecompSVD svd(A);
        Bool_t ok;
        r = svd.Solve(b, ok);
        // if(DEBUG)printf(" ConicalFit : %d %d %lf %lf %lf\n",iter,nus,r[0],r[1],r[3]);
        if (!ok)
            break;
        if ((r[0] * r[0] + r[1] * r[1]) > 1)
        {
            // printf(" ConicalFit >1: iter=%d Nhit=%d %lf %lf %lf\n",iter,nus,r[0],r[1],r[3]);
            break;
        }
        if (ok)
        {
            Flag = 1;
            this->NfiltED = nus;
            sigma = sqrt(sigma / (nus - 3.));
            // sigma = sqrt(sigma);
            dirl = r[0];
            dirm = r[1];
            dirn = sqrt(1 - (r[0] * r[0] + r[1] * r[1]));
            a_temp = AA * alpha + (1 - AA) * r[2];
            // printf("C0 = %lf\n", r[3]);
        }
    }
    // printf("Over\n");
    // get the direction
    phi_temp = -10.;
    the = -10.;
    dt = r[0] * r[0] + r[1] * r[1];
    if (Flag > 0.5 && dt <= 1. && dt >= 0.)
    {
        // printf("sigma = %lf\n", sigma);
        the = asin(sqrt(dt));
        phi_temp = atan2(r[1], r[0]);
        if (phi_temp < 0.)
            phi_temp += DPI * 2.;
        dir[0] = dirl;
        dir[1] = dirm;
        dir[2] = dirn;
        if (AA < 0.5)
            this->a = r[2];
        else
            this->a = alpha;
        this->c0 = r[3];
        this->theta = the;
        this->phi = phi_temp;
        this->Sigma = sigma;
    }
    if (the > -1)
        return 0;
    return -1;
}

double BSTimeoffCalibration::GetMPValue(int id, double delta)
{
    TH1D *h1 = BSTimeoffCalibration::h2->ProjectionX("h1", id, id);
    double mean = h1->GetMean();
    double temp = 0;
    int ite = 30;
    for (int i = 0; i < ite; i++)
    {
        temp = mean;
        h1->SetAxisRange(mean - delta, mean + delta);
        mean = h1->GetMean();
        if (fabs(mean - temp) < 0.01)
            break;
    }
    return mean;
}

int BSTimeoffCalibration::planarfit(TClonesArray &tHits, int np, const char *style)
{
    if (np < 4)
        return -1;
    int i, j, k, nus, Flag = 0;
    double dt, w, c[3], rr[3] = {0, 0, 0};
    double x, y, z, the, phi_temp, sigma, sigmasum, costheta = 0.94, widt = 10.;
    TMatrixD A(3, 3);
    TVectorD b(3), r(3);
    LHHit *tHit;
    for (int iter = 0; iter < 15; iter++)
    {
        A.Zero();
        nus = 0;
        sigmasum = 0;
        for (i = 0; i < np; i++)
        {
            tHit = (LHHit *)((tHits)[i]);
            if (tHit->GetStatus() < 2.5 || tHit->GetPe() < 1.e-3)
                continue; // not use the noise hits
            if (geom->Getxyz(tHit->GetId(), x, y, z, 1, style) < 0)
                continue;
            if (iter == 0)
                w = 1.;
            else
            {
                dt = tHit->GetTime() - ToffsetED[tHit->GetId()] - (rr[0] * x + rr[1] * y + rr[2] - costheta * z) / 0.2998; // ns
                w = exp(-0.5 * pow(dt / widt, 2));                                                                         // sig=10. is better than 5,15,20
                // printf("id = %d, dt = %lf, w = %lf   ", tHit->GetId(), dt, w);
                sigmasum += dt * dt;
                if (iter > 8)
                {
                    if (fabs(dt) > 3 * sigma && fabs(dt) > 30)
                        tHit->SetStatus(1);
                }
            }
            nus++;
            // w*=sqrt(tHit->GetNp()); not good
            // w=sqrt(w); there is a very little change for this weight
            c[0] = x * w;
            c[1] = y * w;
            c[2] = 1 * w;
            for (k = 0; k < 3; k++)
            {
                for (j = 0; j < 3; j++)
                    A[j][k] += c[j] * c[k];
                b[k] += ((tHit->GetTime() - ToffsetED[tHit->GetId()]) * 0.2998 + costheta * z) * w * c[k];
            }
        }
        // solve the equation
        if (nus < 3)
            printf("%d\n", nus);
        if (nus < 3)
            break; // the minimu required hit number is 4
        // printf("%d\n",nus);
        TDecompSVD svd(A);
        Bool_t ok;
        r = svd.Solve(b, ok);
        // if (DEBUG)printf(" PlanarFit : %d %d %lf %lf %lf\n", iter, nus, r[0], r[1], r[2]);
        if (!ok)
        {
            // printf("Can't Solve\n");
            break;
        }
        if ((r[0] * r[0] + r[1] * r[1]) > 1)
        {
            // printf(" PlanarFit >1: %d %d %lf %lf %lf\n",iter,nus,r[0],r[1],r[2]);
            break;
        }
        if (ok)
        {
            sigma = sqrt(sigmasum / (nus - 2.));
            Flag = 1;
            costheta = sqrt(1 - (r[0] * r[0] + r[1] * r[1]));
            rr[0] = r[0];
            rr[1] = r[1];
            rr[2] = r[2];
        }

        // printf("l = %lf, m = %lf, c0 = %lf\n", rr[0], rr[1], rr[2]);
    }
    // get the direction
    phi_temp = -10.;
    the = -10.;
    dt = rr[0] * rr[0] + rr[1] * rr[1];
    if (Flag > 0.5 && dt <= 1. && dt >= 0.)
    {
        the = asin(sqrt(dt));
        phi_temp = atan2(rr[1], rr[0]);
        if (phi_temp < 0.)
            phi_temp += DPI * 2.;
        this->theta = the;
        this->phi = phi_temp;
        this->dir[0] = rr[0];
        this->dir[1] = rr[1];
        this->dir[2] = cos(the);
        // if (DEBUG)printf("%lf %lf\n", the * 57.3, phi_temp * 57.3);
    }

    if (the > -1)
        return 0;
    else
        return -1;
}

double BSTimeoffCalibration::GetRefTime(int IsSet)
{
    int i, id, nus = 0;
    LHHit *tHit;
    double x, y, z, dx, dy, dz, dr, dt;
    double dirl, dirm, dirn;
    double a_temp, c0_temp = 0;
    double reftime, dt_min = 1000;

    TClonesArray *HitsE;
    int ned = lhevent->GetNhitE();
    HitsE = lhevent->GetHitsE();

    if (IsSet == 0)
    {
        dirl = dir[0];
        dirm = dir[1];
        dirn = dir[2];
        a_temp = this->a;
        c0_temp = this->c0;
    }
    else
    {
        dirl = dir_set[0];
        dirm = dir_set[1];
        dirn = dir_set[2];
        a_temp = this->a;
        c0_temp = this->c0_set;
    }

    for (i = 0; i < ned; i++)
    {
        tHit = (LHHit *)((*HitsE)[i]);
        if (tHit->GetStatus() < 2.5 || tHit->GetPe() < 1.e-3)
            continue;
        if (geom->Getxyz(tHit->GetId(), x, y, z, 1, "ED") < 0)
            continue;
        dx = x - core[0];
        dy = y - core[1];
        dz = z;
        dr = (dx * dx + dy * dy + dz * dz - pow(dx * dirl + dy * dirm - dz * dirn, 2.));
        if (dr < 0)
            dr = 0;
        dr = sqrt(dr);

        // dt = tHit->GetTime() - ToffsetED[tHit->GetId()] - (dir[0] * x + dir[1] * y - dir[2] * z + c0) / 0.2998;
        dt = tHit->GetTime() - ToffsetED[tHit->GetId()] - (dirl * x + dirm * y - dirn * z + a_temp * dr + c0_temp) / 0.2998;

        if (nus == 0)
        {
            reftime = dt;
            dt_min = fabs(dt);
        }
        if (fabs(dt) < dt_min)
        {
            dt_min = fabs(dt);
            reftime = dt;
        }
        nus++;
    }
    // printf("Reftime = %lf\n", reftime);
    return reftime;
}

int BSTimeoffCalibration::ApplyPreToffset()
{
    int i, id;
    LHHit *tHit;

    TClonesArray *HitsE;
    int ned = lhevent->GetNhitE();
    HitsE = lhevent->GetHitsE();

    for (i = 0; i < ned; i++)
    {
        tHit = (LHHit *)((*HitsE)[i]);
        if (tHit->GetStatus() < 0 || tHit->GetPe() < 1.e-3)
            continue;
        id = tHit->GetId();
        if (Pre_ToffsetED[id] < -500)
        {
            tHit->SetStatus(-10);
            continue;
        }
        tHit->SetTime(tHit->GetTime() - Pre_ToffsetED[id]);
    }

    return 0;
}

int BSTimeoffCalibration::GetNfilter()
{
    TClonesArray *HitsE;
    int ned = lhevent->GetNhitE();
    HitsE = lhevent->GetHitsE();
    LHHit *tHit;
    int i, nus = 0;

    for (i = 0; i < ned; i++)
    {
        tHit = (LHHit *)((*HitsE)[i]);
        if (tHit->GetStatus() == 5)
            nus++;
    }

    return nus;
}

int BSTimeoffCalibration::eventcallineED_dir()
{
    TClonesArray *HitsE;
    int ned = lhevent->GetNhitE();
    HitsE = lhevent->GetHitsE();

    SetBadHit(*HitsE, ned, "ED");
    Setnparticle(*HitsE, ned, 40.);

    int tWind, rWind;
    int NfiltE;
    tWind = 400;
    rWind = 100;

    NfiltE = spacetimefilter(*HitsE, ned, tWind, rWind, "ED");

    if (!(NfiltE >= 3 && planarfit(*HitsE, ned, "ED") == 0))
        return -1;

    if (conicalfit_dir(*HitsE, ned, 0.035, "ED") != 0)
        return -1;

    NfiltE = noisefilter(*HitsE, ned, -50, 100, rWind + 100, "ED");

    if (NfiltE > 200)
    {
        noisefilter(*HitsE, ned, -50, 100, rWind + 300, "ED");
        if (conicalfit_dir(*HitsE, ned, -0.035, "ED") != 0)
            return -1;
    }
    else
    {
        if (conicalfit_dir(*HitsE, ned, 0.035, "ED") != 0)
            return -1;
    }

    return 0;
}

int BSTimeoffCalibration::ApplyDaz()
{
    Dazcos[0] /= Nevents;
    Dazcos[1] /= Nevents;
    Dazcos[2] /= Nevents;

    printf("dl = %lf, dm = %lf, dn = %lf\n", Dazcos[0], Dazcos[1], Dazcos[2]);

    double x, y, z;
    for (int i = 0; i < MNED; i++)
    {
        if (cal_statED[i] == 0)
            continue;
        if (geom->Getxyz(i, x, y, z, 1, "ED") < 0)
            continue;
        ToffsetED[i] += (Dazcos[0] * x + Dazcos[1] * y - Dazcos[2] * z) / 0.2998;
    }

    return 0;
}

int BSTimeoffCalibration::eventcallineED()
{
    TClonesArray *HitsE;
    int ned = lhevent->GetNhitE();
    HitsE = lhevent->GetHitsE();

    SetBadHit(*HitsE, ned, "ED");

    Setnparticle(*HitsE, ned, 40.);

    int tWind, rWind;
    int NfiltE;
    tWind = 400;
    rWind = 100;

    spacetimefilter(*HitsE, ned, tWind, rWind, "ED");
    if (conicalfit(*HitsE, ned, 0.035, "ED") != 0)
        return -1;
    NfiltE = noisefilter(*HitsE, ned, -50, 100, rWind + 100, "ED");
    if (NfiltE > 200)
    {
        noisefilter(*HitsE, ned, -50, 100, rWind + 300, "ED");
        if (conicalfit(*HitsE, ned, -0.035, "ED") != 0)
            return -2;
    }
    else
    {
        if (conicalfit(*HitsE, ned, 0.035, "ED") != 0)
            return -2;
    }

    Addresidual(*HitsE, ned, "ED");
    // GenerateToffset("ED");

    return 0;
}
