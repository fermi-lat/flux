//  class rootEnergyHist
//  Author:  Theodore Hierath
//  This class contains a histogram and graphing class for
//  angles.
//  The bins are numbered from 0 to number of bins - 1;



#include "rootAngleHist.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>


rootAngleHist::rootAngleHist(int bins) : num_bins(bins)
{
    theta_hist = new rootHist(bins);
    theta_sigma_hist = new rootHist(bins);
    theta_raw_hist = new rootHist(bins);
    phi_hist = new rootHist(bins);
    phi_sigma_hist = new rootHist(bins);
    phi_raw_hist = new rootHist(bins);
    currentType = linear;
    graphTitle = "No Title";
    PhiXLabel = "";
    PhiYLabel = "";
    ThetaXLabel = "";
    ThetaYLabel = "";
}

rootAngleHist::~rootAngleHist(void)
{
    delete theta_hist;
    delete theta_sigma_hist;
    delete theta_raw_hist;
    delete phi_hist;
    delete phi_sigma_hist;
    delete phi_raw_hist;
}

rootAngleHist::rootAngleHist(const rootAngleHist& oldHist) : num_bins(oldHist.num_bins)
{
    theta_hist = new rootHist(num_bins);
    theta_sigma_hist = new rootHist(num_bins);
    theta_raw_hist = new rootHist(num_bins);
    phi_hist = new rootHist(num_bins);
    phi_sigma_hist = new rootHist(num_bins);
    phi_raw_hist = new rootHist(num_bins);
    
    for(int i = 0; i < num_bins; i++)
    {
        double oldCount = (oldHist.theta_hist)->retrieveBin(i);
        theta_hist->updateBin(i, oldCount);
        
        oldCount = (oldHist.theta_sigma_hist)->retrieveBin(i);
        theta_sigma_hist->updateBin(i, oldCount);
        
        oldCount = (oldHist.theta_raw_hist)->retrieveBin(i);
        theta_raw_hist->updateBin(i, oldCount);
        
        oldCount = (oldHist.phi_hist)->retrieveBin(i);
        phi_hist->updateBin(i, oldCount);
        
        oldCount = (oldHist.phi_sigma_hist)->retrieveBin(i);
        phi_sigma_hist->updateBin(i, oldCount);
        
        oldCount = (oldHist.phi_raw_hist)->retrieveBin(i);
        phi_raw_hist->updateBin(i, oldCount);
    }
    currentType = oldHist.currentType;
    graphTitle = oldHist.graphTitle;
    ThetaXLabel = oldHist.ThetaXLabel;
    ThetaYLabel = oldHist.ThetaYLabel;
    PhiXLabel = oldHist.PhiXLabel;
    PhiYLabel = oldHist.PhiYLabel;
}

void rootAngleHist::setTitle(std::string title)
{
    graphTitle = title;
}

void rootAngleHist::setPhiXLabel(std::string label)
{
    PhiXLabel = label;
}

void rootAngleHist::setPhiYLabel(std::string label)
{
    PhiYLabel = label;
}

void rootAngleHist::setThetaXLabel(std::string label)
{
    ThetaXLabel = label;
}

void rootAngleHist::setThetaYLabel(std::string label)
{
    ThetaYLabel = label;
}

void rootAngleHist::setGraphType(const char *graph_type)
{
    if(0 == strcmp(graph_type,"linear") || 0 == strcmp(graph_type,"semilogx"))
        currentType = linear;
    else if(0 == strcmp(graph_type,"log") || 0 == strcmp(graph_type,"semilogy"))
        currentType = log;
    else
        std::cerr << "ERROR:  Invalid Graph Type" << std::endl;
} 

void rootAngleHist::storeTheta(double cos_theta_value)
{
    cos_theta_value += 1;
    
    
    int binNumber = floor((cos_theta_value)/2 * num_bins);
     theta_raw_hist->incrementBin(binNumber);
}

void rootAngleHist::storePhi(double phi_value)
{
    int binNumber = floor(phi_value / (2*M_PI) * num_bins);
    phi_raw_hist->incrementBin(binNumber);
}

double rootAngleHist::retrievePhiCount(int binNumber)
{
    return phi_raw_hist->retrieveBin(binNumber);
}

double rootAngleHist::retrievePhiFlux(int binNumber)
{
    return phi_hist->retrieveBin(binNumber);
}

double rootAngleHist::retrievePhiFluxUncertainty(int binNumber)
{
    return phi_sigma_hist->retrieveBin(binNumber);
}

double rootAngleHist::retrieveThetaCount(int binNumber)
{
    return theta_raw_hist->retrieveBin(binNumber);
}

double rootAngleHist::retrieveThetaFlux(int binNumber)
{
    return theta_hist->retrieveBin(binNumber);
}

double rootAngleHist::retrieveThetaFluxUncertainty(int binNumber)
{
    return theta_sigma_hist->retrieveBin(binNumber);
}

double rootAngleHist::retrieveThetaAngle(int binNumber)
{
    return (binNumber+0.5) / num_bins * 2 - 1;
}

double rootAngleHist::retrievePhiAngle(int binNumber)
{
    return ((binNumber+0.5) / num_bins * 2 * M_PI);
}

// This function must be called before the drawing function
void rootAngleHist::apply(double scale_factor)
{
    for(int i = 0; i < num_bins; ++i)
    {
        double current_flux = theta_hist->retrieveBin(i);
        double current_count = theta_raw_hist->retrieveBin(i);
        current_flux += (current_count * scale_factor);
        theta_hist->updateBin(i,current_flux);
        
        // calculate uncertainty by addition under quadrature
        current_flux = sqrt(pow(theta_sigma_hist->retrieveBin(i),2)+current_count) * scale_factor;
        theta_sigma_hist->updateBin(i,current_flux);
        
        // clear temporary array
        theta_raw_hist->updateBin(i,0);
        
        current_flux = phi_hist->retrieveBin(i);
        current_count = phi_raw_hist->retrieveBin(i);
        current_flux += (current_count * scale_factor);
        phi_hist->updateBin(i,current_flux);
        
        // calculate uncertainty
        current_flux = sqrt(pow(phi_sigma_hist->retrieveBin(i),2)+current_count) * scale_factor;
        phi_sigma_hist->updateBin(i,current_flux);
        
        // clear temporary array
        phi_raw_hist->updateBin(i,0);
    }
}

void rootAngleHist::reset(void)
{
    for(int i = 0; i < num_bins; ++i)
    {
        phi_hist->updateBin(i,0);
        phi_sigma_hist->updateBin(i,0);
        phi_raw_hist->updateBin(i,0);
        theta_hist->updateBin(i,0);
        theta_sigma_hist->updateBin(i,0);
        theta_raw_hist->updateBin(i,0);
    }
    
    graphTitle = "No Title";
    PhiXLabel = "";
    PhiYLabel = "";
    ThetaXLabel = "";
    ThetaYLabel = "";
}

void rootAngleHist::draw(double scale_factor, std::string mode, int current_plot, int total_plots)
{
    char *window_title = "Graph Window";
    
    if(current_plot >= total_plots)
    {
        std::cerr << "Error:  Invalid plot number" << std::endl;
        return;
//      exit(0);
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
            "   char *angle_window_title = \"" << window_title << "\";\n"
            "   char *theta_graph_title = \"" << "Particle Flux vs. Zenith Angle" << "\";\n"
            "   char *phi_graph_title = \"" << "Particle Flux vs. Azimuth Angle" << "\";\n"
            "   char *theta_y_label = \"" << ThetaYLabel << "\";\n"
            "   char *theta_x_label = \"" << ThetaXLabel << "\";\n"
            "   char *phi_y_label = \"" << PhiYLabel << "\";\n"
            "   char *phi_x_label = \"" << PhiXLabel << "\";\n"
            "   int num_bins = "<< num_bins << ";\n"
            "   double theta_angle[num_bins];\n"
            "   double e_theta_angle[num_bins];\n"
            "   double phi_angle[num_bins];\n"
            "   double e_phi_angle[num_bins];\n"
            
            "   for(int i = 0; i < num_bins; i++) {\n"
            "      theta_angle[i] = (i+0.5) / num_bins * 2 - 1;\n"
            "      e_theta_angle[i] = 0;\n"
            "      phi_angle[i] = (double(i)+0.5)/double(num_bins)*360;\n"
            "      e_phi_angle[i] = 0;\n"
            "   }\n"
            
            "   //                 name, title, wtopx, wtopy, ww, wh\n"
            "   c2 = new TCanvas(\"c2\",angle_window_title, 200, 200, 700, 500);\n"
            "   c2->SetGrid();\n"
            "   c2->GetFrame()->SetFillColor(21);\n"
            "   c2->GetFrame()->SetBorderSize(12);\n"
            "   c3 = new TCanvas(\"c3\",angle_window_title, 200, 200, 700, 500);\n"
            "   c3->SetGrid();\n"
            "   c3->GetFrame()->SetFillColor(21);\n"
            "   c3->GetFrame()->SetBorderSize(12);\n";
        
        out_file << 
            "   theta_leg = new TLegend(0.73,0.83,0.99,0.99);\n"
            "   phi_leg = new TLegend(0.73,0.83,0.99,0.99);\n"
            "   double scale_factor" << current_plot << " = " << scale_factor << ";\n"
            "   double theta_count" << current_plot << "[] = {\n";
        {for(int i = 0; i < num_bins; i++)
            out_file << std::setw(5) << retrieveThetaFlux(i) << (i%5==4? ",\n" : ",");
        }
        
        out_file <<
            "   };\n"
            "   double phi_count" << current_plot << "[] = {\n";
        {for(int i = 0; i < num_bins; i++)
            out_file << std::setw(5) << retrievePhiFlux(i) << (i%5==4? ",\n" : ",");
        }
        
        out_file << 
            "   };\n"
            "   double upper_angle;"
            "   double lower_angle;"
            "   double e_theta_count" << current_plot << "[] = {\n";
        {for(int i = 0; i < num_bins; i++)
            out_file << std::setw(5) << retrieveThetaFluxUncertainty(i) << (i%5==4? ",\n" : ",");
        }
        out_file <<
            "  };\n"
            "   double e_phi_count" << current_plot << "[] = {\n";
        {for(int i = 0; i < num_bins; i++)
            out_file << std::setw(5) << retrievePhiFluxUncertainty(i) << (i%5==4? ",\n" : ",");
        }
        out_file <<
            "   };\n"
            "   for(int i = 0; i < num_bins; i++) {\n"
            "      e_theta_count" << current_plot << "[i] *= scale_factor" << current_plot << ";\n"
            "      e_phi_count" << current_plot << "[i] *= scale_factor" << current_plot << ";\n"
            "      theta_count" << current_plot << "[i] *= scale_factor" << current_plot << ";\n"
            "      phi_count" << current_plot << "[i] *= scale_factor" << current_plot << ";\n"
            "   }\n"
            "   c2->cd();"
            "   theta_gr" << current_plot << " = new TGraphErrors(num_bins,theta_angle,theta_count" << current_plot 
            << ",e_theta_angle,e_theta_count" << current_plot << ");\n"
            "   theta_gr" << current_plot << "->SetTitle(theta_graph_title);\n"
            "   theta_gr" << current_plot << "->SetMarkerColor(2);\n"
            "   theta_gr" << current_plot << "->SetMarkerStyle(21);\n"
            "   theta_leg->AddEntry(theta_gr" << current_plot << ",\"" << graphTitle << "\",\"P\");\n"
            "   c3->cd();"
            "   phi_gr" << current_plot << " = new TGraphErrors(num_bins,phi_angle,phi_count" << current_plot
            << ",e_phi_angle,e_phi_count" << current_plot << ");\n"
            "   phi_gr" << current_plot << "->SetTitle(phi_graph_title);\n"
            "   phi_gr" << current_plot << "->SetMarkerColor(2);\n"
            "   phi_gr" << current_plot << "->SetMarkerStyle(21);\n"
            "   phi_leg->AddEntry(phi_gr" << current_plot << ",\"" << graphTitle << "\",\"P\");\n";
        
        out_file.close();
   }   
   else if(current_plot <= total_plots - 1)
   {
       std::ofstream out_file("graph.cxx", std::ios::app);
       out_file <<
           "   double scale_factor" << current_plot << " = " << scale_factor << ";\n"
           "   double theta_count" << current_plot << "[] = {\n";
       {for(int i = 0; i < num_bins; i++)
           out_file << std::setw(5) << retrieveThetaFlux(i) << (i%5==4? ",\n" : ",");
       }
       
       out_file <<
           "   };\n"
           "   double phi_count" << current_plot << "[] = {\n";
       {for(int i = 0; i < num_bins; i++)
           out_file << std::setw(5) << retrievePhiFlux(i) << (i%5==4? ",\n" : ",");
       }
       
       out_file << 
           "   };\n"
           "   double e_theta_count" << current_plot << "[] = {\n";
       {for(int i = 0; i < num_bins; i++)
           out_file << std::setw(5) << retrieveThetaFluxUncertainty(i) << (i%5==4? ",\n" : ",");
       }
       out_file <<
           "   };\n"
           "   double e_phi_count" << current_plot << "[] = {\n";
       {for(int i = 0; i < num_bins; i++)
           out_file << std::setw(5) << retrievePhiFluxUncertainty(i) << (i%5==4? ",\n" : ",");
       }
       out_file <<
           "   };\n"
           "   for(i = 0; i < num_bins; i++) {\n"
           "      e_theta_count" << current_plot << "[i] *= scale_factor" << current_plot << ";\n"
           "      e_phi_count" << current_plot << "[i] *= scale_factor" << current_plot << ";\n"
           "      theta_count" << current_plot << "[i] *= scale_factor" << current_plot << ";\n"
           "      phi_count" << current_plot << "[i] *= scale_factor" << current_plot << ";\n"
           "   }\n"
           "   c2->cd();"
           "   theta_gr" << current_plot << " = new TGraphErrors(num_bins,theta_angle,theta_count" << current_plot 
           << ",e_theta_angle,e_theta_count" << current_plot << ");\n"
           "   theta_gr" << current_plot << "->SetMarkerColor(" << 2+current_plot << ");\n"
           "   theta_gr" << current_plot << "->SetMarkerStyle(" << 21 << ");\n"
           "   theta_leg->AddEntry(theta_gr" << current_plot << ",\"" << graphTitle << "\",\"P\");\n"
           "   c3->cd();"
           "   phi_gr" << current_plot << " = new TGraphErrors(num_bins,phi_angle,phi_count" << current_plot
           << ",e_phi_angle,e_phi_count" << current_plot << ");\n"
           "   phi_gr" << current_plot << "->SetMarkerColor(" << 2+current_plot << ");\n"
           "   phi_gr" << current_plot << "->SetMarkerStyle(" << 21 << ");\n"
           "   phi_leg->AddEntry(phi_gr" << current_plot << ",\"" << graphTitle << "\",\"P\");\n";
       out_file.close();
   }
   
   if(current_plot >= total_plots - 1)
   {
       std::ofstream out_file("graph.cxx", std::ios::app);
       
       out_file << 
           "   c2->cd();"
           "   double theta_max_count = 0;\n"
           "   double theta_min_count = 1e12;\n"
           "   for(int i = 0; i < num_bins; i++)\n"
           "   {\n";
       {for(int plot = 0; plot < total_plots; plot++) {
           out_file << 
               "      if(theta_count" << plot << "[i] > theta_max_count)\n"
               "         theta_max_count = theta_count" << plot << "[i];\n"
               "      if(theta_count" << plot << "[i] < theta_min_count && theta_count" << plot << "[i] > 0)\n"
               "         theta_min_count = theta_count" << plot << "[i];\n";
       }}
       out_file <<
           "   }\n"
           "   double theta_angle_limits[] = {-1,1};\n"
           "   double theta_count_limits[] = {0,theta_max_count};\n"
           "   int theta_bin_limits = 2;\n"
           "   theta_graph0 = new TGraph(theta_bin_limits,theta_angle_limits,theta_count_limits);\n"
           "   theta_graph0->SetTitle(theta_graph_title);\n"
           "   theta_graph0->Draw(\"AP\");\n"
           "   TAxis *theta_ax = theta_graph0->GetXaxis();\n"
           "   TAxis *theta_ay = theta_graph0->GetYaxis();\n"
           "   theta_ay->SetTitle(theta_y_label); theta_ay->CenterTitle(1);\n"
           "   theta_ax->SetLimits(-1,1);\n"
           "   theta_ax->SetTitle(theta_x_label); theta_ax->CenterTitle(1); \n"
           "   theta_ax->SetTitleOffset(1.2);\n";
       
       {for(int plot = 0; plot < total_plots; plot++) {
           out_file <<
               "   theta_gr" << plot << "->Draw(\"P\");\n";
       }}
       
       out_file << 
           "   theta_leg->Draw();\n"
           "   c2->Modified();\n"
           "   c2->Update();\n";
       
       out_file << 
           "   c3->cd();"
           "   double phi_max_count = 0;\n"
           "   double phi_min_count = 1e12;\n"
           "   for(int i = 0; i < num_bins; i++)\n"
           "   {\n";
       {for(int plot = 0; plot < total_plots; plot++) {
           out_file << 
               "      if(phi_count" << plot << "[i] > phi_max_count)\n"
               "         phi_max_count = phi_count" << plot << "[i];\n"
               "      if(phi_count" << plot << "[i] < phi_min_count && phi_count" << plot << "[i] > 0)\n"
               "         phi_min_count = phi_count" << plot << "[i];\n";
       }}
       out_file <<
           "   }\n"
           "   double phi_angle_limits[] = {0,360};\n"
           "   double phi_count_limits[] = {0,phi_max_count};\n"
           "   int phi_bin_limits = 2;\n"
           "   phi_graph0 = new TGraph(phi_bin_limits,phi_angle_limits,phi_count_limits);\n"
           "   phi_graph0->SetTitle(phi_graph_title);\n"
           "   phi_graph0->Draw(\"AP\");\n"
           "   TAxis *phi_ax = phi_graph0->GetXaxis();\n"
           "   TAxis *phi_ay = phi_graph0->GetYaxis();\n"
           "   phi_ay->SetTitle(phi_y_label); phi_ay->CenterTitle(1);\n"
           "   phi_ax->SetLimits(0,360);\n"
           "   phi_ax->SetTitle(phi_x_label); phi_ax->CenterTitle(1); \n"
           "   phi_ax->SetTitleOffset(1.2);\n";
       
       {for(int plot = 0; plot < total_plots; plot++) {
           out_file <<
               "   phi_gr" << plot << "->Draw(\"P\");\n";
       }}
       
       out_file << 
           "   phi_leg->Draw();\n"
           "   c3->Modified();\n"
           "   c3->Update();\n";
       
       if(mode != "begin")
           out_file << "\n}\n";
       
       out_file.close();
       
       if(mode != "begin")
           system("root -l graph.cxx");
   }
}

void rootAngleHist::writeFile(double scale_factor, std::ostream& out_file)
{
    
}
