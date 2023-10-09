// CREATED BY LUCIAN SCHARENBERG

void fellowStyle(int axisStyle = 0, double defaultPadWidth = 1000, double defaultPadHeight = 1000) {

    // Different axes styles
    // 0 -> default:       left and bottom
    // 1 -> second y axis: left, right and bottom
    // 1 -> 2D:            left, bottom and COLZ
        
    // Default white background for all plots
    gROOT->SetStyle("Plain");
   
    double stylePadWidth;
    double stylePadHeight;

    double defaultMarginLeft = 0.20;
    double defaultMarginBottom = 0.20;
    double defaultMarginRight = 0.07;
    double defaultMarginTop = 0.07;

    double marginLeft;
    double marginRight;
    double marginTop;
    double marginBottom;

    double defaultTitleOffsetX = 1.85;
    double defaultTitleOffsetY = 1.85;

    double titleOffsetX;
    double titleOffsetY;

    if (axisStyle == 0) {

	stylePadWidth = defaultPadWidth;
	stylePadHeight = defaultPadHeight;
	marginLeft = defaultMarginLeft;
	marginRight = defaultMarginRight;
	marginTop = defaultMarginTop;
	marginBottom = defaultMarginBottom;
	titleOffsetX = defaultTitleOffsetX;
	titleOffsetY = defaultTitleOffsetY;

    }

    else if (axisStyle == 1) {

	double colzFactor = 1.2;
	stylePadWidth = colzFactor * defaultPadWidth;
	stylePadHeight = defaultPadHeight;
	marginLeft = defaultMarginLeft * defaultPadWidth / stylePadWidth;
	marginRight = (defaultMarginRight + (stylePadWidth - defaultPadWidth) / defaultPadWidth) * defaultPadWidth / stylePadWidth;
	marginTop = defaultMarginTop;
	marginBottom = defaultMarginBottom;
	titleOffsetX = defaultTitleOffsetX;
	titleOffsetY = defaultTitleOffsetY * defaultPadWidth / stylePadWidth;

    }

    else if (axisStyle == 2) {

	double padFactor = 1.36;
	stylePadWidth = defaultPadWidth;
	stylePadHeight = padFactor * defaultPadHeight;
	marginLeft = defaultMarginLeft * defaultPadWidth / stylePadWidth;
	marginRight = defaultMarginRight;
	marginTop = defaultMarginTop;
	marginBottom = defaultMarginBottom;
	titleOffsetX = defaultTitleOffsetX * defaultPadWidth / stylePadWidth;
	titleOffsetY = defaultTitleOffsetY;

    }
    
    else if (axisStyle == 3) {

	double axisFactor = 1.13; // 0.2 new margin - 0.07 original margin
	stylePadWidth = axisFactor * defaultPadWidth;
	stylePadHeight = defaultPadHeight;
	marginLeft = defaultMarginLeft * defaultPadWidth / stylePadWidth;
	marginRight = (defaultMarginRight + (stylePadWidth - defaultPadWidth) / defaultPadWidth) * defaultPadWidth / stylePadWidth;
	marginTop = defaultMarginTop;
	marginBottom = defaultMarginBottom;
	titleOffsetX = defaultTitleOffsetX;
	titleOffsetY = defaultTitleOffsetY * defaultPadWidth / stylePadWidth;

    }

    gStyle->SetCanvasDefW(stylePadWidth);
    gStyle->SetCanvasDefH(stylePadHeight);

    gStyle->SetPadBottomMargin(marginBottom);
    gStyle->SetPadTopMargin(marginTop);
    gStyle->SetPadRightMargin(marginRight);
    gStyle->SetPadLeftMargin(marginLeft);

    gStyle->SetCanvasColor(0);
    gStyle->SetStatColor(0);
    gStyle->SetPadColor(0);
    gStyle->SetTitleFillColor(0);

    // Label and title offset
    gStyle->SetTitleOffset(titleOffsetX, "x");
    gStyle->SetTitleOffset(titleOffsetY, "yz");

    double labelOffset = 0.015;
    gStyle->SetLabelOffset(labelOffset, "x"); // * titleOffsetX / defaultTitleOffsetX, "x");
    gStyle->SetLabelOffset(labelOffset, "yz"); // * titleOffsetY / defaultTitleOffsetY, "yz");

    // Use Helvetica as text font
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "xyz");
    gStyle->SetTitleFont(42, "xyz");
    gStyle->SetStatFont(42);
    gStyle->SetLegendFont(42);

    gStyle->SetLabelColor(kBlack, "xyz");

    // Set the text size of the title and labels
    double textSize = 0.055;
    gStyle->SetLabelSize(textSize, "xyz");
    gStyle->SetTitleSize(textSize, "xyz");

    //gStyle->SetLineWidth(2);
    gStyle->SetPadTickX(1);
    if (axisStyle == 3) {
        gStyle->SetPadTickY(0);
    }
    else {
	gStyle->SetPadTickY(1);
    }
    gStyle->SetFrameLineWidth(1);
    gStyle->SetFrameLineColor(kBlack);
    gStyle->SetAxisColor(kBlack, "xyz");

    // Set the number of divisions to show
    gStyle->SetNdivisions(505, "xyz");
    
    // Set ticksize
    double tickLengthX = 0.035 * stylePadWidth / defaultPadWidth;
    double tickLengthY = 0.035 * defaultPadWidth / stylePadWidth;
    gStyle->SetTickLength(tickLengthX, "x");
    gStyle->SetTickLength(tickLengthY, "yz");

    // Histograms
    gStyle->SetDrawBorder(0);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetHistLineWidth(2);

    // Fits
    gStyle->SetOptFit(0);
    gStyle->SetFuncWidth(2);

    // Data points
    gStyle->SetEndErrorSize(0);
    gStyle->SetMarkerStyle(20);
    gStyle->SetMarkerSize(2);

    // Other stuff
    gStyle->SetPalette(kBird);
    gStyle->SetNumberContours(500);
    TGaxis::SetMaxDigits(4);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetGridColor(18);
    gStyle->SetGridStyle(1);
    
    gROOT->ForceStyle();
   
}
