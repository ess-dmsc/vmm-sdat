// Obtain the spatial resolution and other parameters
// --------------------------------------------------
// lucian.scharenberg@cern.ch
// July and August 2023

using namespace ROOT;
using namespace std;


double getEfficiency(const char *filename, const char *histHit, const char *histNoHit) {

    TFile *f = new TFile(filename);
    TH2D *hHit = (TH2D*)f->Get(histHit);
    TH2D *hNoHit = (TH2D*)f->Get(histNoHit);

    double nHit = hHit->GetEntries();
    double nNoHit = hNoHit->GetEntries();

    return (nHit / (nHit + nNoHit));
    
}

double getEfficiencyError(const char *filename, const char *histHit, const char *histNoHit) {

    TFile *f = new TFile(filename);
    TH2D *hHit = (TH2D*)f->Get(histHit);
    TH2D *hNoHit = (TH2D*)f->Get(histNoHit);

    double nHit = hHit->GetEntries();
    double nNoHit = hNoHit->GetEntries();
    
    double a = nHit;
    double b = nNoHit;
    double aErr = TMath::Sqrt(a);
    double bErr = TMath::Sqrt(b);

    double dda = b * aErr / ((a + b) * (a + b));
    double ddb = - a * bErr / ((a + b) * (a + b));

    return TMath::Sqrt(dda * dda + ddb * ddb);
    
}

TH1D *getHistogramAnamicom(const char *filename, const char *histname) {

    TFile *f = new TFile(filename);
    TH1D *h = (TH1D*)f->Get(histname);
    h->SetDirectory(0);

    return h;

}

pair<double, double> getAmplitude(const char *filename, const char *histname) {

    TFile *f = new TFile(filename);
    TH1D *h = (TH1D*)f->Get(histname);
    h->SetDirectory(0);

    return {h->GetMean(), h->GetMeanError()};
  
}

pair<double, double> getClusterSize(const char *filename, const char *histname) {

    TFile *f = new TFile(filename);
    TH1D *h = (TH1D*)f->Get(histname);
    h->SetDirectory(0);

    return {h->GetMean(), h->GetMeanError()};
  
}

double fitDoubleGaussian(double *x, double *par) {

    double mean, sigmaCore, sigmaTail, a, scale;
    double coreGauss, tailGauss;

    mean = par[0];
    sigmaCore = par[1];
    sigmaTail = par[2];
    a = par[3];
    scale = par[4];

    coreGauss = a * TMath::Exp(- TMath::Power(x[0] - mean, 2) / (2 * TMath::Power(sigmaCore, 2)));
    tailGauss = (1 - a) * TMath::Exp(- TMath::Power(x[0] - mean, 2) / (2 * TMath::Power(sigmaTail, 2)));

    return (scale * (coreGauss + tailGauss));

}

double getSigma(double a, double sigmaCore, double sigmaTail) {

    return TMath::Sqrt(a * TMath::Power(sigmaCore, 2) + (1 - a) * TMath::Power(sigmaTail, 2));

}

double getSigmaError(double a, double aError, double sigmaCore,
                     double sigmaCoreError, double sigmaTail, double sigmaTailError) {

    double aDerivative = (TMath::Power(sigmaCore, 2) - TMath::Power(sigmaTail, 2)) / (2 * TMath::Sqrt(a * (TMath::Power(sigmaCore, 2) - TMath::Power(sigmaTail, 2)) + TMath::Power(sigmaTail, 2)));
    double sigmaCoreDerivative = (a * sigmaCore) / (2 * TMath::Sqrt(a * (TMath::Power(sigmaCore, 2) - TMath::Power(sigmaTail, 2)) + TMath::Power(sigmaTail, 2)));
    double sigmaTailDerivative = (- (a - 1) * sigmaTail) / (2 * TMath::Sqrt(a * (TMath::Power(sigmaCore, 2) - TMath::Power(sigmaTail, 2)) + TMath::Power(sigmaTail, 2)));

    double aPart = TMath::Power(aDerivative * aError, 2);
    double sigmaCorePart = TMath::Power(sigmaCoreDerivative * sigmaCoreError, 2);
    double sigmaTailPart = TMath::Power(sigmaTailDerivative * sigmaTailError, 2);

    return TMath::Sqrt(aPart + sigmaCorePart + sigmaTailPart);

}

double geometricMean(double sigmaIn, double sigmaOut) {

    return (TMath::Sqrt(sigmaIn * sigmaOut));

}

double geometricMeanError(double sigmaIn, double sigmaInError, double sigmaOut, double sigmaOutError) {

    double a = sigmaOutError * sigmaIn / (2 * geometricMean(sigmaIn, sigmaOut));
    double b = sigmaInError * sigmaOut / (2 * geometricMean(sigmaIn, sigmaOut));

    return TMath::Sqrt(TMath::Power(a, 2) + TMath::Power(b, 2));

}

double trackResolution(double z, double zT1, double zT2, double zT3, double resT1, double resT2, double resT3) {

    vector<double> position = {zT1, zT2, zT3};
    vector<double> weights = {resT1, resT2, resT3};

    double l11 = 0.0, l12 = 0.0, l22 = 0.0;
    double a, b;

    for (int i = 0; i < position.size(); i++) {
        l11 += 1.0 / TMath::Power(weights[i], 2);
        l12 += position[i] / TMath::Power(weights[i], 2);
        l22 += TMath::Power(position[i], 2) / TMath::Power(weights[i], 2);
    }

    a = l22 - 2 * z * l12 + z * z * l11;
    b = l11 * l22 - l12 * l12;

    return TMath::Sqrt(a / b);

}

double trackResolutionError(double z, double zT1, double zT2, double zT3, double resT1, double resT2, double resT3,
                            double zError, double zT1Error, double zT2Error, double zT3Error, double resT1Error, double resT2Error, double resT3Error) {

    vector<double> position = {zT1, zT2, zT3};
    vector<double> weights = {resT1, resT2, resT3};
    vector<double> weightsError = {resT1Error, resT2Error, resT3Error};

    double l11 = 0.0, l12 = 0.0, l22 = 0.0;
    double l11Error = 0.0, l12Error = 0.0, l22Error = 0.0;
    double l11Derivative, l12Derivative, l22Derivative;
    double zDerivative = 0.0;

    for (int i = 0; i < position.size(); i++) {
        l11 += 1.0 / TMath::Power(weights[i], 2);
        l11Error += TMath::Power((-2.0 * weightsError[i] / TMath::Power(weights[i], 3)), 2);
        l12 += position[i] / TMath::Power(weights[i], 2);
        l12Error += (TMath::Power(-2.0 * z * weightsError[i] / TMath::Power(weights[i], 3), 2) +
                     TMath::Power(zError / (weights[i] * weights[i]), 2));
        l22 += TMath::Power(position[i], 2) / TMath::Power(weights[i], 2);
        l22Error += (TMath::Power(-2.0 * z * z * weightsError[i] / TMath::Power(weights[i], 3), 2) +
                     TMath::Power(2 * z * zError / (weights[i] * weights[i]), 2));
    }

    l11Error = TMath::Sqrt(l11Error);
    l12Error = TMath::Sqrt(l12Error);
    l22Error = TMath::Sqrt(l22Error);

    double a = l22;
    double b = l12;
    double c = l11;

    double sqrtPlaceholder = TMath::Sqrt((a - 2 * b * z + c * z * z) / (a * c - b * b));

    l22Derivative = ((1 / (a * c - b * b) - (c * (a - 2 * b * z + c * z * z)) / TMath::Power(a * c - b * b, 2))) / (2 * sqrtPlaceholder);
    l12Derivative = (2 * b * (a - 2 * b * z + c * z * z) / TMath::Power(a * c - b * b, 2) - (2 * z) / (a * c - b * b)) / (2 * sqrtPlaceholder);
    l11Derivative = ((z * z / (a * c - b * b)) - (a * (a - 2 * b * z + c * z * z)) / TMath::Power(a * c - b * b, 2)) / (2 * sqrtPlaceholder);
    zDerivative = (2 * c * z - 2 * b) / (2 * (a * c - b * b) * sqrtPlaceholder);

    return TMath::Sqrt(TMath::Power(l22Derivative * l22Error, 2) +
                       TMath::Power(l12Derivative * l12Error, 2) +
                       TMath::Power(l11Derivative * l11Error, 2) +
                       TMath::Power(zDerivative * zError, 2));

}

double getResolutionDUT(double sigmaDetector, double sigmaTrack) {

    return TMath::Sqrt(sigmaDetector * sigmaDetector - sigmaTrack * sigmaTrack);

}

double getResolutionDUTError(double sigmaDetector, double sigmaDetectorError,
                             double sigmaTrack, double sigmaTrackError) {

    double a = sigmaDetector * sigmaDetectorError / getResolutionDUT(sigmaDetector, sigmaTrack);
    double b = -1 * sigmaTrack * sigmaTrackError / getResolutionDUT(sigmaDetector, sigmaTrack);

    return TMath::Sqrt(a * a + b * b);

}

double myGaus(double *x, double *par) {

    double scale, mean, sigma;

    mean = par[1];
    sigma = par[2];
    scale = par[0];

    return scale * TMath::Exp(- TMath::Power(x[0] - mean, 2) / (2 * TMath::Power(sigma, 2)));

}

pair<double, double> getWidth(TH1D *h, string name, string runNumber, int nBins, string rms) {
  
    auto c = new TCanvas("c", "", 800, 600);

    TH1D *originalPlot = new TH1D(*h);
    TH1D *plot = dynamic_cast<TH1D*>(originalPlot->Rebin(nBins, "plot"));
    double sigma = 9999.0;
    double sigmaError = 9999.0;

    if (rms == "fullFit") {

	plot->GetXaxis()->SetRangeUser(-0.6, 0.6);
	plot->GetXaxis()->SetTitle("Residuals / mm");
	plot->GetYaxis()->SetTitle("Counts");
    
	TF1 *fit = new TF1("fit", fitDoubleGaussian, -1.0, 1.0, 5);
	fit->SetParameters(plot->GetMean(), plot->GetRMS() * 0.5, plot->GetRMS() * 1.5, 0.8, plot->GetBinContent(plot->GetMaximumBin()));
	fit->SetParNames("mean", "sigmaCore", "sigmaTail", "a", "scale");
	fit->SetParLimits(1, 0.0, 1.0);
	fit->SetParLimits(2, 0.0, 10.0);
	fit->SetParLimits(3, 0.5, 1.0);
	fit->SetLineColor(kRed);
	plot->Fit("fit", "BL", "", -0.6, 0.6);

	double mean = fit->GetParameter("mean");
	double meanError = fit->GetParError(0);
	double sigmaCore = fit->GetParameter("sigmaCore");
	double sigmaCoreError = fit->GetParError(1);
	double sigmaTail = fit->GetParameter("sigmaTail");
	double sigmaTailError = fit->GetParError(2);
	double a = fit->GetParameter("a");
	double aError = fit->GetParError(3);
	double scale = fit->GetParameter("scale");
	double scaleError = fit->GetParError(4);
	
	sigma = getSigma(a, sigmaCore, sigmaTail);
	sigmaError = getSigmaError(a, aError, sigmaCore, sigmaCoreError, sigmaTail, sigmaTailError);
	//cout << "SIGMA: "<< sigma << " +- " << sigmaError << endl;

	TF1 *fitCore = new TF1("fitCore", myGaus, -1.0, 1.0, 3);
	fitCore->SetParameters(scale * a, mean, sigmaCore);
	fitCore->SetLineColor(kRed);
	fitCore->SetLineStyle(kDashed);
	fitCore->Draw("SAME");

	TF1 *fitTail = new TF1("fitTail", myGaus, -1.0, 1.0, 3);
	fitTail->SetParameters(scale * (1.0 - a), mean, sigmaTail);
	fitTail->SetLineColor(kRed);
	fitTail->SetLineStyle(kDashed);
	fitTail->Draw("SAME");

    }
    else if (rms == "simpleFit") {
	
	plot->GetXaxis()->SetRangeUser(-0.6, 0.6);
	plot->GetXaxis()->SetTitle("Residuals / mm");
	plot->GetYaxis()->SetTitle("Counts");
    
	TF1 *fit = new TF1("fit", "gaus");
	fit->SetParameters(plot->GetBinContent(plot->GetMaximumBin()), plot->GetMean(), plot->GetRMS() * 0.5);
	fit->SetParNames("scale", "mean", "sigma");
	fit->SetParLimits(2, 0.0, 1.0);
	fit->SetLineColor(kRed);

	string inName = "resinfit0";
	string outName = "resoutfit0";
	if (name.find(outName) != string::npos) {
	    plot->Fit("fit", "BL", "", -0.07, 0.07);
	}
	else if (name.find(inName) != string::npos) {
	    plot->Fit("fit", "BL", "", -0.03, 0.03);
	}

	sigma = fit->GetParameter("sigma");
	sigmaError = fit->GetParError(2);
	
    }
    else if (rms == "rms") {

	plot->GetXaxis()->SetRangeUser(-1.0, 1.0);
	plot->GetXaxis()->SetTitle("Residuals / mm");
	plot->GetYaxis()->SetTitle("Counts");

	sigma = plot->GetRMS();
	sigmaError = plot->GetRMSError();
	cout << sigma << endl;

    }

    c->Update();
    c->SaveAs(("./checkFit/run" + runNumber + "_" + name + "_rebin_" + to_string(nBins) + ".pdf").c_str());
    delete c;
    return {sigma, sigmaError};

}

pair<double, double> getTrackResolution(string fileName, string runNumber,
                                        string indexT1, string indexT2, string indexT3,
                                        double zT1, double zT2, double zT3,
                                        double zDUT, double zError,
                                        int nBins, string rms) {

    // Get the tracking detector residual widths (excluded from fit)
    pair<double, double> widthOutT1 = getWidth(getHistogramAnamicom(fileName.c_str(), ("resoutfit0mm" + indexT1).c_str()), ("resoutfit0mm" + indexT1).c_str(), runNumber, nBins, rms);
    pair<double, double> widthOutT2 = getWidth(getHistogramAnamicom(fileName.c_str(), ("resoutfit0mm" + indexT2).c_str()), ("resoutfit0mm" + indexT2).c_str(), runNumber, nBins, rms);
    pair<double, double> widthOutT3 = getWidth(getHistogramAnamicom(fileName.c_str(), ("resoutfit0mm" + indexT3).c_str()), ("resoutfit0mm" + indexT3).c_str(), runNumber, nBins, rms);

    // Get the tracking detector residual widths (included in fit)
    pair<double, double> widthInT1 = getWidth(getHistogramAnamicom(fileName.c_str(), ("resinfit0mm" + indexT1).c_str()), ("resinfit0mm" + indexT1).c_str(), runNumber, nBins, rms);
    pair<double, double> widthInT2 = getWidth(getHistogramAnamicom(fileName.c_str(), ("resinfit0mm" + indexT2).c_str()), ("resinfit0mm" + indexT2).c_str(), runNumber, nBins, rms);
    pair<double, double> widthInT3 = getWidth(getHistogramAnamicom(fileName.c_str(), ("resinfit0mm" + indexT3).c_str()), ("resinfit0mm" + indexT3).c_str(), runNumber, nBins, rms);

    // Calculate the actual resolution
    double resT1 = geometricMean(widthInT1.first, widthOutT1.first);
    double resT2 = geometricMean(widthInT2.first, widthOutT2.first);
    double resT3 = geometricMean(widthInT3.first, widthOutT3.first);

    // Calculate the resolution error
    double resT1Error = geometricMeanError(widthInT1.first, widthInT1.second, widthOutT1.first, widthOutT1.second);
    double resT2Error = geometricMeanError(widthInT2.first, widthInT2.second, widthOutT2.first, widthOutT2.second);
    double resT3Error = geometricMeanError(widthInT3.first, widthInT3.second, widthOutT3.first, widthOutT3.second);

    // Calculate the track resolution
    double result = trackResolution(zDUT, zT1, zT2, zT3, resT1, resT2, resT3);
    double resultError = trackResolutionError(zDUT, zT1, zT2, zT3, resT1, resT2, resT3, zError, zError, zError, zError, resT1Error, resT2Error, resT3Error);

    return {result, resultError};

}

void spatialResolution() {

    // Decide which method to use to determine the width of the
    // residual distributions
    // * "fullFit" (default)
    //   -> fit a double Gaussian function, as described in my thesis
    // * "simpleFit"
    //   -> fit only the central part
    // * "rms"
    //   -> take the RMS of the histogram (no fit at all, depends
    //      strongly on the axis range, i.e. not recommended)
    string rms = "simpleFit";

    // Decide which method to use to determine the spatial resolution
    // * "full" (default)
    //   -> fit everything and subtract the track uncertainty, as
    //      described in my thesis
    // * "geometric"
    //   -> just use the geometric mean as described in here
    //      https://doi.org/10.1016/j.nima.2004.08.132
    string resMethod = "full";

    // Decide which parameter, in addition to the spatial resolution
    // is supposed to be plotted.
    // All values of
    // * Detector plane efficiency        -> 'E'
    // * Average cluster charge per plane -> 'A'
    // * Average cluster size per plane   -> 'S'
    // are calculated and stored in the CSV file, but only one of them
    // is plotted. The default is the efficiency.
    char plotParameter = 'E';
  
    // Detector positions 
    double zT1 = 0.0;
    double zT2 = -550.0;
    double zT3 = -960.0;
    double zDUT = -359.5;
    double zError = 2.0;

    // Detector indices
    string iT1X = "0";
    string iT2X = "3";
    string iT3X = "5";

    string iT1Y = "1";
    string iT2Y = "2";
    string iT3Y = "4";

    string iDUTX = "6";
    string iDUTY = "7";

    // In case rebinning of the residual distribution is needed,
    // change 1 to the value you would like to have
    int nBins = 1;    
    
    // Data vectors
    vector<double> voltage = {30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 600.0, 900.0, 1100.0, 1300.0, 1500.0, 1700.0};
    vector<double> gain;
    vector<double> gainError;

    vector<double> resolutionX;
    vector<double> resolutionY;
    vector<double> resolutionXError;
    vector<double> resolutionYError;

    vector<double> efficiencyX;
    vector<double> efficiencyY;
    vector<double> efficiencyXError;
    vector<double> efficiencyYError;

    vector<double> amplitudeX;
    vector<double> amplitudeY;
    vector<double> amplitudeXError;
    vector<double> amplitudeYError;

    vector<double> clusterSizeX;
    vector<double> clusterSizeY;
    vector<double> clusterSizeXError;
    vector<double> clusterSizeYError;
    
    for (int g = 0; g < voltage.size(); g++) {
        gain.push_back(voltage[g] / 0.3);
        gainError.push_back(0.0);
    }
    
    // Files to loop over
    string scanName = "tb202307";
    string plotNameBase = "softwareTest";
    string plotName = (string)("./" + plotNameBase + nBins + ".pdf");
    vector<string> runNumbers = {"54", "55", "56", "57", "58", "59", "60", "61", "62", "63", "64", "65", "66", "67", "68", "69"};

    vector<string> fileList;
    for (string item : runNumbers) {
        if (item.size() == 1) {
	    item = "000" + item;
        }
        else if (item.size() == 2) {
	    item = "00" + item;
        }
        else if (item.size() == 3) {
            item = "0" + item;
        }
        fileList.push_back("/eos/project/r/rd51/SRSVMMData/" + scanName + "/trackdata/anarun" + item + ".root");
    }

    // Get the resolution and the efficiency
    double res, resError, eff, effError, ampl, amplError, size, sizeError;

    for (int i = 0; i < runNumbers.size(); i++) {

	if (resMethod == "full") {

	    pair<double, double> trackingResults;
	    pair<double, double> detectorResults;

	    trackingResults = getTrackResolution(fileList[i].c_str(), runNumbers[i], iT1X, iT2X, iT3X, zT1, zT2, zT3, zDUT, zError, nBins, rms);
	    detectorResults = getWidth(getHistogramAnamicom(fileList[i].c_str(), ("resoutfit0mm" + iDUTX).c_str()), ("resoutfit0mm" + iDUTX).c_str(), runNumbers[i], nBins, rms);
	    res = getResolutionDUT(detectorResults.first, trackingResults.first);
	    resError = getResolutionDUTError(detectorResults.first, detectorResults.second, trackingResults.first, trackingResults.second);
	    resolutionX.push_back(res * 1000.0);
	    resolutionXError.push_back(resError * 1000.0);
    
	    trackingResults = getTrackResolution(fileList[i].c_str(), runNumbers[i], iT1Y, iT2Y, iT3Y, zT1, zT2, zT3, zDUT, zError, nBins, rms);
	    detectorResults = getWidth(getHistogramAnamicom(fileList[i].c_str(), ("resoutfit0mm" + iDUTY).c_str()), ("resoutfit0mm" + iDUTY).c_str(), runNumbers[i], nBins, rms);
	    res = getResolutionDUT(detectorResults.first, trackingResults.first);
	    resError = getResolutionDUTError(detectorResults.first, detectorResults.second, trackingResults.first, trackingResults.second);
	    resolutionY.push_back(res * 1000.0);
	    resolutionYError.push_back(resError * 1000.0);

	}
	else if (resMethod == "geometric") {

	    pair<double, double> detectorResultsOut;
	    pair<double, double> detectorResultsIn;
	    
	    detectorResultsOut = getWidth(getHistogramAnamicom(fileList[i].c_str(), ("resoutfit0mm" + iDUTX).c_str()), ("resoutfit0mm" + iDUTX).c_str(), runNumbers[i], nBins, rms);
	    detectorResultsIn = getWidth(getHistogramAnamicom(fileList[i].c_str(), ("resinfit0mm" + iDUTX).c_str()), ("resinfit0mm" + iDUTX).c_str(), runNumbers[i], nBins, rms);
	    res = geometricMean(detectorResultsIn.first, detectorResultsOut.first);
	    resError = geometricMeanError(detectorResultsIn.first, detectorResultsIn.second, detectorResultsOut.first, detectorResultsOut.second);
	    resolutionX.push_back(res * 1000.0);
	    resolutionXError.push_back(resError * 1000.0);

	    detectorResultsOut = getWidth(getHistogramAnamicom(fileList[i].c_str(), ("resoutfit0mm" + iDUTY).c_str()), ("resoutfit0mm" + iDUTY).c_str(), runNumbers[i], nBins, rms);
	    detectorResultsIn = getWidth(getHistogramAnamicom(fileList[i].c_str(), ("resinfit0mm" + iDUTY).c_str()), ("resinfit0mm" + iDUTY).c_str(), runNumbers[i], nBins, rms);
	    res = geometricMean(detectorResultsIn.first, detectorResultsOut.first);
	    resError = geometricMeanError(detectorResultsIn.first, detectorResultsIn.second, detectorResultsOut.first, detectorResultsOut.second);
	    resolutionY.push_back(res * 1000.0);
	    resolutionYError.push_back(resError * 1000.0);

	}

	eff = getEfficiency(fileList[i].c_str(), ("hits0mm" + iDUTX).c_str(), ("inefficspots0mm" + iDUTX).c_str());
	effError = getEfficiencyError(fileList[i].c_str(), ("hits0mm" + iDUTX).c_str(), ("inefficspots0mm" + iDUTX).c_str());
	efficiencyX.push_back(eff * 100.0);
	efficiencyXError.push_back(effError * 100.0);

	eff = getEfficiency(fileList[i].c_str(), ("hits0mm" + iDUTY).c_str(), ("inefficspots0mm" + iDUTY).c_str());
	effError = getEfficiencyError(fileList[i].c_str(), ("hits0mm" + iDUTY).c_str(), ("inefficspots0mm" + iDUTY).c_str());
	efficiencyY.push_back(eff * 100.0);
	efficiencyYError.push_back(effError * 100.0);

	ampl = getAmplitude(fileList[i].c_str(), ("totcha" + iDUTX).c_str()).first;
	amplError = getAmplitude(fileList[i].c_str(), ("totcha" + iDUTX).c_str()).second;
	amplitudeX.push_back(ampl);
	amplitudeXError.push_back(amplError);

	ampl = getAmplitude(fileList[i].c_str(), ("totcha" + iDUTY).c_str()).first;
	amplError = getAmplitude(fileList[i].c_str(), ("totcha" + iDUTY).c_str()).second;
	amplitudeY.push_back(ampl);
	amplitudeYError.push_back(amplError);

	size = getClusterSize(fileList[i].c_str(), ("nohitstrips" + iDUTX).c_str()).first;
	sizeError = getClusterSize(fileList[i].c_str(), ("nohitstrips" + iDUTX).c_str()).second;
	clusterSizeX.push_back(size);
	clusterSizeXError.push_back(sizeError);

	size = getClusterSize(fileList[i].c_str(), ("nohitstrips" + iDUTY).c_str()).first;
	sizeError = getClusterSize(fileList[i].c_str(), ("nohitstrips" + iDUTY).c_str()).second;
	clusterSizeY.push_back(size);
	clusterSizeYError.push_back(sizeError);
        
    }

    // Save the data to file
    ofstream csvFile("./dataPoints/" + plotNameBase + nBins + ".csv");
    if (csvFile.is_open()) {
	csvFile << "Gain"                    << "," << "Gain Error" << ","
		<< "Resolution X (um)"       << "," << "Resolution X Error (um)"       << ","
		<< "Resolution Y (um)"       << "," << "Resolution Y Error (um)"       << ","
		<< "Efficiency X (%)"        << "," << "Efficiency X Error (%)"        << ","
		<< "Efficiency Y (%)"        << "," << "Efficiency Y Error (%)"        << ","
		<< "Mean charge X (ADC)"     << "," << "Mean charge X Error (ADC)"     << ","
		<< "Mean charge Y (ADC)"     << "," << "Mean charge Y Error (ADC)"     << ","
		<< "Cluster size X (Strips)" << "," << "Cluster size X Error (Strips)" << ","
		<< "Cluster size Y (Strips)" << "," << "Cluster size Y Error (Strips)" << endl;
        for (int line = 0; line < resolutionX.size(); line++) {
            csvFile << gain[line] << ","  << gainError[line] << ","
                    << resolutionX[line]  << "," << resolutionXError[line]  << ","
                    << resolutionY[line]  << "," << resolutionYError[line]  << ","
                    << efficiencyX[line]  << "," << efficiencyXError[line]  << ","
                    << efficiencyY[line]  << "," << efficiencyYError[line]  << ","
		    << amplitudeX[line]   << "," << amplitudeXError[line]   << ","
                    << amplitudeY[line]   << "," << amplitudeYError[line]   << ","
		    << clusterSizeX[line] << "," << clusterSizeXError[line] << ","
                    << clusterSizeY[line] << "," << clusterSizeYError[line] << endl;
        }
    }
    else {
        cout << "CSV file could not be opened" << endl;
    }
    
    // Plot the data
    auto cFull = new TCanvas("cFull", "", 800, 600);
    cFull->cd();
  
    auto mgRes = new TMultiGraph();

    double limitX1 = 0.0;
    double limitX2 = 6000.0;
    double limitY1 = 0.0;
    double limitY2 = 140.0;
    double limitSecondY1 = 0.0;
    double limitSecondY2 = 0.0;

    switch (plotParameter) {

    case 'E':

	limitSecondY1 = 0.0;
	limitSecondY2 = 100.0;
	break;
	
    case 'A':

	limitSecondY1 = 0.0;
	limitSecondY2 = 1500.0;
	break;

    case 'S':

	limitSecondY1 = 0.0;
	limitSecondY2 = 10.0;
	break;
	
    }
    
    auto graphX = new TGraphErrors((int)resolutionX.size(), &gain[0], &resolutionX[0], &gainError[0], &resolutionXError[0]);
    graphX->SetLineColor(kBlack);
    graphX->SetMarkerColor(kBlack);
    graphX->SetMarkerStyle(20);
    mgRes->Add(graphX);

    auto graphY = new TGraphErrors((int)resolutionY.size(), &gain[0], &resolutionY[0], &gainError[0], &resolutionYError[0]);
    graphY->SetLineColor(kRed);
    graphY->SetMarkerColor(kRed);
    graphY->SetMarkerStyle(21);
    mgRes->Add(graphY);

    auto mgOther = new TMultiGraph();

    // Scale the effAmpl to match the second Y axis
    double scalingFactor = (limitY2 - limitY1) / (limitSecondY2 - limitSecondY1) + limitSecondY1;

    for (int i = 0; i < gain.size(); i++) {

	switch (plotParameter) {
	    
	case 'E':
	    
	    efficiencyX[i] = efficiencyX[i] * scalingFactor;
	    efficiencyY[i] = efficiencyY[i] * scalingFactor;
	    efficiencyXError[i] = efficiencyXError[i] * scalingFactor;
	    efficiencyYError[i] = efficiencyYError[i] * scalingFactor;
	    break;

	case 'A':

	    amplitudeX[i] = amplitudeX[i] * scalingFactor;
	    amplitudeY[i] = amplitudeY[i] * scalingFactor;
	    amplitudeXError[i] = amplitudeXError[i] * scalingFactor;
	    amplitudeYError[i] = amplitudeYError[i] * scalingFactor;
	    break;

	case 'S':

	    clusterSizeX[i] = clusterSizeX[i] * scalingFactor;
	    clusterSizeY[i] = clusterSizeY[i] * scalingFactor;
	    clusterSizeXError[i] = clusterSizeXError[i] * scalingFactor;
	    clusterSizeYError[i] = clusterSizeYError[i] * scalingFactor;
	    break;

	}

    }

    TGraphErrors *graphOtherX;
    TGraphErrors *graphOtherY;
    
    switch (plotParameter) {

    case 'E':
	    
	graphOtherX = new TGraphErrors((int)efficiencyX.size(), &gain[0], &efficiencyX[0], &gainError[0], &efficiencyXError[0]);
	graphOtherY = new TGraphErrors((int)efficiencyY.size(), &gain[0], &efficiencyY[0], &gainError[0], &efficiencyYError[0]);
	break;
    
    case 'A':

	graphOtherX = new TGraphErrors((int)amplitudeX.size(), &gain[0], &amplitudeX[0], &gainError[0], &amplitudeXError[0]);
	graphOtherY = new TGraphErrors((int)amplitudeY.size(), &gain[0], &amplitudeY[0], &gainError[0], &amplitudeYError[0]);
	break;

    case 'S':

	graphOtherX = new TGraphErrors((int)clusterSizeX.size(), &gain[0], &clusterSizeX[0], &gainError[0], &clusterSizeXError[0]);
	graphOtherY = new TGraphErrors((int)clusterSizeY.size(), &gain[0], &clusterSizeY[0], &gainError[0], &clusterSizeYError[0]);
	break;
	    
    }

    graphOtherX->SetLineColor(kBlack);
    graphOtherX->SetMarkerColor(kBlack);
    graphOtherX->SetMarkerStyle(24);
    mgOther->Add(graphOtherX);

    graphOtherY->SetLineColor(kRed);
    graphOtherY->SetMarkerColor(kRed);
    graphOtherY->SetMarkerStyle(25);
    mgOther->Add(graphOtherY);

    mgRes->GetXaxis()->SetTitle("Drift field / V/cm");
    mgRes->GetXaxis()->SetLimits(limitX1, limitX2);
    mgRes->GetYaxis()->SetTitle("Spatial resolution / #mum");
    mgRes->GetYaxis()->SetRangeUser(limitY1, limitY2);
    mgRes->Draw("APL");
    mgOther->Draw("PL");
  
    TGaxis *axis = new TGaxis(limitX2, limitY1, limitX2, limitY2, limitSecondY1, limitSecondY2, 510, "+L");
    axis->SetLabelFont(42);
    axis->SetTitleFont(42);

    auto legend = new TLegend(0.6, 0.45, 0.85, 0.65);
    legend->SetTextFont(42);
    legend->SetBorderSize(0.0);
    legend->AddEntry(graphX, "Res. X", "lep");
    legend->AddEntry(graphY, "Res. Y", "lep");

    switch (plotParameter) {
    
    case 'E':

	axis->SetTitle("Efficiency / %");
	legend->AddEntry(graphOtherX, "Eff. X", "lep");
	legend->AddEntry(graphOtherY, "Eff. Y", "lep");
	break;

    case 'A':
	
	axis->SetTitle("Mean signal charge / ADC counts");
	legend->AddEntry(graphOtherX, "Ampl. X", "lep");
	legend->AddEntry(graphOtherY, "Ampl. Y", "lep");
	break;

    case 'S':

	axis->SetTitle("Average cluster size / Strips");
	legend->AddEntry(graphOtherX, "Size X", "lep");
	legend->AddEntry(graphOtherY, "Size Y", "lep");
	break;
	    
    }
    
    axis->Draw();
    legend->Draw();
  
    cFull->Update();
    cFull->SaveAs("./spatialResolution.pdf");
    //cFull->SaveAs(plotName.c_str());
    
}
