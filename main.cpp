#include <iostream>
#include <TH1D.h>
#include <TTree.h>
#include <TFile.h>
#include "string"
#include "TCanvas.h"
#include <fstream>
#include <sstream>
#include <numeric>
#include <TF1.h>
#include <TStyle.h>

std::vector<std::string> getNextLineAndSplitIntoTokens(std::istream &str){
    /// This method reads in a csv line by line and separates everything
    std::vector<std::string>    result;
    std::string            line;
    std::getline(str, line);

    std::stringstream linestream(line);
    std::string            cell;

    while(std::getline(linestream, cell, ' '))
        if(!cell.empty()) result.push_back(cell);

    if(!linestream && cell.empty()) result.emplace_back("");

    return result;
}

int main() {
    //This program aims to correct for the faulty calibration of the XYTER-Chip by integration

    const int num_com = 31;  // Number of comparators
    const std::vector<std::array<std::string, 3>> settings{
        {"../data/outputcal_201119_143946.txt",             // calibration location
         "../data/XYTER_201119_143946_for7200s.root",       // 241 Am spectrum location
         "../results/241Am_201119_143946.pdf"},             // output file location
        {"../data/outputcal_201120_101026.txt",             // calibration location
         "../data/XYTER_201120_101026_for60s.root",       // 241 Am spectrum location
         "../results/241Am_201120_101026.pdf"}};            // output file location

    const int s_no = 1;
    //Get the weights of each comparator
    std::ifstream myfile;
    myfile.open (settings.at(s_no).at(0));
    assert(myfile.is_open());

    for(int i=0; i<32; i++) auto temp = getNextLineAndSplitIntoTokens(myfile);

    std::vector<std::vector<int>> raw_calibration;  // first dimension comparator (1 to 31), second: charges
    for(int i=0; i<num_com; i++) raw_calibration.emplace_back();

    for(std::vector<std::string> temp = getNextLineAndSplitIntoTokens(myfile); temp.size() != 1;
        temp = getNextLineAndSplitIntoTokens(myfile)){
        for(int i=0; i<num_com; i++){
            raw_calibration.at(i).push_back(std::stoi(temp.at(32-i)));
        }
    }

    // After having all S-curves in binary form, we can simply build the sum to integrate
    std::vector<int> sum_comp;
    for(auto &i: raw_calibration) sum_comp.push_back(std::accumulate(i.begin(), i.end(), 0));

    std::vector<double> weights(num_com + 1);
    weights.at(0) = 0.5;
    weights.back() = num_com + 0.5;
    for (int i=1; i<31; i++) weights.at(i) = std::max(0.1, (double)
            (sum_comp.at(i-1) - sum_comp.at(i))/(sum_comp.at(0)-sum_comp.back()) *
            (num_com - 1)) + weights.at(i-1);

    //weights.at(15) = 15.8;

    // Read in the root file containing the 241Am source data
    auto spectrumfile = TFile::Open(settings.at(s_no).at(1).c_str());
    assert(spectrumfile);
    auto spectrumtree = (TTree*)spectrumfile->Get("otree");
    assert(spectrumtree);

    // Get a distribution of events in a nice std::vector
    int ADC;
    spectrumtree->SetBranchAddress("nadc",&ADC);
    TH1D raw_hist("raw241AmSource", "{}^{241}Am Source test, 2h", num_com, 0.5, num_com + 0.5);
    int nentries = (Int_t)spectrumtree->GetEntries();
    for (int i=0; i<nentries; i++) {
        spectrumtree->GetEntry(i);
        raw_hist.Fill(ADC);
    }
    std::vector<int> raw_spectrum;
    for(int i = 1; i < (num_com + 1); i++) raw_spectrum.push_back(raw_hist.GetBinContent(i));

    // Merge data into transformed histogram.
    TH1D cleaned_hist("241AmSource", "{}^{241}Am Source test, 2h", num_com, weights.data());
    for(int i = 1; i < (num_com + 1); i++){
        double bin_width = weights.at(i)-weights.at(i-1);
        if(!bin_width) bin_width += 1;
        std::cout << "Bin Width for ADC " << i << " = " << bin_width << std::endl;
        cleaned_hist.SetBinContent(i, (double)raw_hist.GetBinContent(i)/bin_width);
    }

    TF1 *f1 = new TF1("f1","gaus",11,30);
    cleaned_hist.Fit(f1, "R");

    TCanvas source_canvas("sourcecanvas", "{}^{241}Am Source test, 2h",
                                 800,450);
    gStyle->SetOptFit(0012);
    gStyle->SetStatX(0.41);
    gStyle->SetStatY(0.8);
    gStyle->SetStatH(0.1);
    gStyle->SetStatW(0.14);
    source_canvas.SetTicks();
    cleaned_hist.GetYaxis()->SetTitle("counts");
    cleaned_hist.GetXaxis()->SetTitle("ADC bit");
    cleaned_hist.Draw();
    source_canvas.SaveAs(settings.at(s_no).at(2).c_str(), ".pdf");

    return 0;
}