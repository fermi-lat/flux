//  class rootEnergyHist
//  Author:  Theodore Hierath
//  This class contains a histogram and graphing class for
//  energies.
//  The bins are numbered from 0 to number of bins - 1;

#include "rootEnergyHist.h"
#include <iostream>
#include <iomanip>
#include <fstream>

rootEnergyHist::rootEnergyHist(int bins, 
                               double min_energy, 
                               double max_energy) :
emin(min_energy),
emax(max_energy),
num_bins(bins)
{
    range = log10(emax/emin);
    currentType = linear;
    flux_hist = new rootHist(bins);
    flux_sigma_hist = new rootHist(bins);
    raw_hist = new rootHist(bins);


    graphTitle = "No Title";
    xlabel = "";
    ylabel = "";
    fluxMode = false;
    use_flux_min = false;
    use_flux_max = false;
}

rootEnergyHist::~rootEnergyHist(void)
{
    delete flux_hist;
    delete flux_sigma_hist;
    delete raw_hist;
}

rootEnergyHist::rootEnergyHist(const rootEnergyHist& oldHist) :
emin(oldHist.emin),
emax(oldHist.emax),
num_bins(oldHist.num_bins)
{
    range = oldHist.range;
    currentType = oldHist.currentType;
    flux_hist = new rootHist(num_bins);
    flux_sigma_hist = new rootHist(num_bins);
    raw_hist = new rootHist(num_bins);
    
    for(int i = 0; i < num_bins; i++)
    {  
        double oldCount = (oldHist.raw_hist)->retrieveBin(i);
        raw_hist->updateBin(i, oldCount );
        
        oldCount = (oldHist.flux_hist)->retrieveBin(i);
        flux_hist->updateBin(i, oldCount );
        
        oldCount = (oldHist.flux_sigma_hist)->retrieveBin(i);
        flux_sigma_hist->updateBin(i, oldCount );
    }
    
    graphTitle = oldHist.graphTitle;
    xlabel = oldHist.xlabel;
    ylabel = oldHist.ylabel;
    fluxMode = oldHist.fluxMode;
    use_flux_min = oldHist.use_flux_min;
    use_flux_max = oldHist.use_flux_max;
}


void rootEnergyHist::setTitle(std::string title)
{
    graphTitle = title;
}

void rootEnergyHist::setXLabel(std::string label)
{
    xlabel = label;
}

void rootEnergyHist::setYLabel(std::string label)
{
    ylabel = label;
}

void rootEnergyHist::setGraphType(const char *graph_type)
{
    if(0 == strcmp(graph_type,"linear"))
        currentType = linear;
    else if(0 == strcmp(graph_type,"semilogx"))
        currentType = semilogx;
    else if(0 == strcmp(graph_type,"semilogy"))
        currentType = semilogy;
    else if(0 == strcmp(graph_type,"log"))
        currentType = loglog;
    else
        std::cerr << "ERROR:  Invalid Graph Type" << std::endl;
}

void rootEnergyHist::store(double energy)
{ 
    if(energy < emax && energy > emin)
    {
        if(currentType == loglog || currentType == semilogx)
        {
            int currentBin = int(num_bins * (log10(energy/emin) / range));
            raw_hist->incrementBin(currentBin);
        }
        else
        {
            int currentBin = int(floor(num_bins * (energy - emin) / (emax - emin)));
            raw_hist->incrementBin(currentBin);
        }
    }
}

double rootEnergyHist::retrieveCount(int binNumber)
{
    return raw_hist->retrieveBin(binNumber);
}

double rootEnergyHist::retrieveFlux(int binNumber)
{
    return flux_hist->retrieveBin(binNumber);
}

double rootEnergyHist::retrieveFluxUncertainty(int binNumber)
{
    return flux_sigma_hist->retrieveBin(binNumber);
}

double rootEnergyHist::retrieveEnergy(int binNumber)
{
    return emin*pow(10.0,((binNumber + 0.5)/num_bins) * range );
}

double rootEnergyHist::retrieveRange(void)
{
    return range;
}


void rootEnergyHist::setFluxMode(void)
{
    fluxMode = true;
}

void rootEnergyHist::setFluxMin(double f)
{
    use_flux_min = true;
    flux_min = f;
}

void rootEnergyHist::setFluxMax(double f)
{
    use_flux_max = true;
    flux_max = f;
}

// This function must be called before the drawing function
void rootEnergyHist::apply(double scale_factor)
{
    for(int i = 0; i < num_bins; ++i)
    {
        double current_flux = flux_hist->retrieveBin(i);
        double current_count = raw_hist->retrieveBin(i);
        current_flux += (current_count * scale_factor);
        flux_hist->updateBin(i,current_flux);
        
        // calculate uncertainty by addition under quadrature
        current_flux = sqrt(pow(flux_sigma_hist->retrieveBin(i),2)+current_count) * scale_factor;
        flux_sigma_hist->updateBin(i,current_flux);
        
        // clear temporary array
        raw_hist->updateBin(i,0);
    }
}

void rootEnergyHist::reset(void)
{
    for(int i = 0; i < num_bins; ++i)
    {
        flux_hist->updateBin(i,0);
        flux_sigma_hist->updateBin(i,0);
        raw_hist->updateBin(i,0);
    }
    
    graphTitle = "No Title";
    xlabel = "";
    ylabel = "";
    fluxMode = false;
    use_flux_min = false;
    use_flux_max = false;
    range = log10(emax/emin);
    currentType = linear;
}


// scale_factor = number to multiply each bin by
// mode = "normal", "begin", or "end"
// current_plot = current_plot number (numbering starts at 0)
// total_plots = the total number of plots to go on the graph
void rootEnergyHist::draw(double scale_factor, std::string mode, int current_plot, int total_plots)
{
    char *window_title = "Graph Window";
    
    if(current_plot >= total_plots)
    {
        std::cerr << "Error:  Invalid plot number" << std::endl;
        return;
//        exit(0);
    }
    
    if(current_plot == 0)
    {
        std::ofstream out_file;
        
        if(mode != "end")
        {
            out_file.open("graph.cxx");
            
            if(false == out_file.is_open())
            {
                std::cerr << "Unable to open temporary file for writing." << std::endl;
                return;
//                exit(0);
            }
            
            out_file << 
                "{\n"
                "   gROOT->Reset();\n";
        }
        else
        {
            out_file.open("graph.cxx", std::ios::app);   
        }
        
        out_file <<         
            "   char *energy_window_title = \"" << window_title << "\";\n"
            "   char *energy_graph_title = \"" << "Particle Flux vs. Kinetic Energy" << "\";\n"
            "   char *energy_y_label = \"" << ylabel << "\";\n"
            "   char *energy_x_label = \"" << xlabel << "\";\n"
            
            "   double energy_min= "<< emin << ", energy_max=" << emax << ";\n"
            "   int num_bins = "<< num_bins << ";\n"
            "   double log_energy_range = " << retrieveRange() << ";\n"
            "   double energy[num_bins];\n"
            "   double e_energy[num_bins];\n"
            
            "   for(int i = 0; i < num_bins; i++) {\n"
            "      energy[i] = energy_min*pow(10,((i+0.5)/num_bins) * log_energy_range);\n"
            "      e_energy[i] = 0;\n"
            "   }\n"
            "   //                 name, title, wtopx, wtopy, ww, wh\n"
            "   c1 = new TCanvas(\"c1\",energy_window_title, 200, 200, 700, 500);\n"
            "   c1->SetGrid();\n"
            "   c1->GetFrame()->SetFillColor(21);\n"
            "   c1->GetFrame()->SetBorderSize(12);\n";
        
        if(currentType == loglog)        out_file << "   c1->SetLogx(1); c1->SetLogy(1);\n";
        else if(currentType == semilogx) out_file << "   c1->SetLogx(1);\n";
        else if(currentType == semilogy) out_file << "   c1->SetLogy(1);\n";
        else if(currentType == linear)   ;
        else
            std::cerr << "Error:  invalid graph type" << std::endl;
        
        out_file << 
            "   leg = new TLegend(0.73,0.83,0.99,0.99);\n"
            "   double scale_factor" << current_plot << " = " << scale_factor << ";\n"
            "   double count" << current_plot << "[] = {\n";
        {for(int i = 0; i < num_bins; ++i)
            out_file << std::setw(7) << retrieveFlux(i) << (i%5==4? ",\n" : ",");
        }
        
        out_file << 
            "   };\n"
            "   double e_count" << current_plot << "[num_bins] = {\n";
        {for(int i = 0; i < num_bins; ++i)
            out_file << std::setw(7) << retrieveFluxUncertainty(i) << (i%5==4? ",\n" : ",");
        }
        
        out_file <<
            "   };\n"
            "   for(i = 0; i < num_bins; i++) {\n"
            "      e_count" << current_plot << "[i] *= scale_factor" << current_plot;
        if(fluxMode == true)
            out_file << "/(energy[i]*1000);\n";
        else
            out_file << ";\n";
        
        out_file <<
            "      count" << current_plot << "[i] *=scale_factor" << current_plot;
        
        if(fluxMode == true)
            out_file << "/(energy[i]*1000);\n";
        else
            out_file << ";\n";
        
        out_file <<
            "   }\n"
            "   gr" << current_plot << " = new TGraphErrors(num_bins,energy,count" << current_plot 
            << ",e_energy,e_count" << current_plot << ");\n"
            "   gr" << current_plot << "->SetTitle(energy_graph_title);\n"
            "   gr" << current_plot << "->SetMarkerColor(2);\n"
            "   gr" << current_plot << "->SetMarkerStyle(21);\n"
            "   leg->AddEntry(gr" << current_plot << ",\"" << graphTitle << "\",\"P\");\n";
        
        out_file.close();
    }   
    else if(current_plot <= total_plots - 1)
    {
        std::ofstream out_file("graph.cxx", std::ios::app);
        out_file <<
            "   c1->cd();\n"
            "   double scale_factor" << current_plot << " = " << scale_factor << ";\n"
            "   double count" << current_plot << "[] = {\n";
        {for(int i = 0; i < num_bins; i++)
            out_file << std::setw(5) << retrieveFlux(i) << (i%5==4? ",\n" : ",");
        }
        
        out_file << 
            "   };\n"
            "   double e_count" << current_plot << "[num_bins] = {\n";
        {for(int i = 0; i < num_bins; ++i)
            out_file << std::setw(7) << retrieveFluxUncertainty(i) << (i%5==4? ",\n" : ",");
        }
        
        out_file <<
            "   };\n"
            "   for(i = 0; i < num_bins; i++) {\n"
            "      e_count" << current_plot << "[i] *= scale_factor" << current_plot;
        if(fluxMode == true)
            out_file << "/(energy[i]*1000);\n";
        else
            out_file << ";\n";
        
        out_file <<
            "      count" << current_plot << "[i] *=scale_factor" << current_plot;
        
        if(fluxMode == true)
            out_file << "/(energy[i]*1000);\n";
        else
            out_file << ";\n";
        
        out_file <<
            "   }\n"
            "   gr" << current_plot << " = new TGraphErrors(num_bins,energy,count" << current_plot 
            << ",e_energy,e_count" << current_plot << ");\n"
            "   gr" << current_plot << "->SetMarkerColor(" << 2+current_plot << ");\n"
            "   gr" << current_plot << "->SetMarkerStyle(" << 21 << ");\n"
            "   leg->AddEntry(gr" << current_plot << ",\"" << graphTitle << "\",\"P\");\n";
        out_file.close();
    }
    
    if(current_plot >= total_plots - 1)
    {
        std::ofstream out_file("graph.cxx", std::ios::app);
        
        out_file << 
            "   double max_count = 0;\n"
            "   double min_count = 1e12;\n"
            "   for(int i = 0; i < num_bins; i++)\n"
            "   {\n";
        {for(int plot = 0; plot < total_plots; plot++) {
            out_file << 
                "      if(count" << plot << "[i] > max_count)\n"
                "         max_count = count" << plot << "[i];\n"
                "      if(count" << plot << "[i] < min_count && count" << plot << "[i] > 0)\n"
                "         min_count = count" << plot << "[i];\n";
        }}
        out_file <<
            "   }\n"
            "   double energy_limits[] = {energy_min,energy_max};\n"
            "   double count_limits[] = {";
        
        if(use_flux_min)
            out_file << flux_min << ",";
        else
            out_file << "min_count,";
        
        if(use_flux_max)
            out_file << flux_max;
        else
            out_file << "max_count";
        
        out_file <<
            "};\n"
            "   int bin_limits = 2;\n"
            "   graph0 = new TGraph(bin_limits,energy_limits,count_limits);\n"
            "   graph0->SetTitle(energy_graph_title);\n"
            "   graph0->Draw(\"AP\");\n"
            "   TAxis *ax = graph0->GetXaxis();\n"
            "   TAxis *ay = graph0->GetYaxis();\n"
            "   ay->SetTitle(energy_y_label); ay->CenterTitle(1);\n"
            "   ax->SetLimits(energy_min, energy_max);\n"
            "   ax->SetTitle(energy_x_label); ax->CenterTitle(1); \n"
            "   ax->SetTitleOffset(1.2);\n";
        
        {for(int plot = 0; plot < total_plots; plot++) {
            out_file <<
                "   gr" << plot << "->Draw(\"P\");\n";
        }}
        
        out_file << 
            "   leg->Draw();\n"
            "   c1->Modified();\n"
            "   c1->Update();\n";
        
        if(mode != "begin")
            out_file << "\n}\n";
        
        out_file.close();
        
        if(mode != "begin")
            system("root -l graph.cxx");
    }
}


// If the output is going to be used for another source, the histogram should have 
// a large enough number of bins to minimize rounding errors caused by bin size.
void rootEnergyHist::writeFile(double scale_factor, std::ostream& out_file)
{
    out_file << "%Energy (in MeV)  vs. Flux (in particles/m2-s-sr-MeV)" << std::endl;
    
    for(int i = 0; i < num_bins; i++)
    {
        // Use center of bin to determine the energies
        double energy = retrieveEnergy(i) * 1000; // GeV->MeV
        double flux = double(retrieveFlux(i)) * scale_factor / energy; // E*Flux -> Flux
        
        out_file << std::setiosflags(std::ios::right|std::ios::scientific) 
            << std::setw(14) << std::setprecision(5) << energy 
            << std::setiosflags(std::ios::right|std::ios::scientific)
            << std::setw(14) << std::setprecision(5) << flux 
            << std::endl;
    }
    
}
