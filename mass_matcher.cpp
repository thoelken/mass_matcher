#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <cmath>

using namespace std;


struct ion {
  double mass;
  string text;
  ion() {};
  ion(double m, string t): mass(m), text(t) {}
};


struct xpsm {
  ion scan;
  map<double, ion> candidates;
};


map<double, ion> parse(const char* filename) {
  
  ifstream file;
  file.open(filename);
  map<double, ion> candidates;
  
  if(file.is_open()) {
  
    string line;
    while ( getline (file,line) ) {
    
      istringstream ls (line);
      string text;
      double mass;
      
      if(!(ls >> mass)) {
        cerr << "could not parse mass value from line";
      }
      while(ls) {
        ls >> text;
      }
      
      candidates.insert(pair<double, ion>(mass, ion(mass, text)));
    }
    file.close();
    
  } else {
    cout << "Unable to open file";
  }
  
  return candidates;
}


map<double, ion> parseMGF(const char* filename) {

	const string PEPMASS = "PEPMASS=";
	const string TITLE = "TITLE=";
	const string CHARGE = "CHARGE=";
	const double PROTON = 1.007276466812;
  
  ifstream file;
  file.open(filename);
  map<double, ion> scans;
  
  if(file.is_open()) {
  
  	string text;
    double mass;
    int charge;
    
    string line;
    while(getline(file, line)) {
    	if(line.length() > PEPMASS.length() && line.compare(0, PEPMASS.length(), PEPMASS) == 0) {
    		istringstream ls (line.substr(PEPMASS.length()));
    		ls >> mass;
    		mass = mass*charge-charge*PROTON;
//    		cout << text << "\t" << mass << "\t" << charge << endl;
    		scans.insert(pair<double, ion>(mass, ion(mass, text)));
    	} else if(line.length() > TITLE.length() && line.compare(0, TITLE.length(), TITLE) == 0) {
    		text = line.substr(TITLE.length(), line.find_last_not_of(" \t\r\n")-TITLE.length()+1);
    	} else if(line.length() > CHARGE.length() && line.compare(0, CHARGE.length(), CHARGE) == 0) {
    		istringstream ls (line.substr(CHARGE.length()));
    		ls >> charge;
    	}
    
    }
    file.close();
    
  } else {
    cout << "Unable to open file";
  }
  
  return scans;
}


xpsm matchScan(map<double, ion> xlinks, ion scan, double error = 0.00001, int top = 100) {

  map<double, ion> submap ( xlinks.lower_bound(scan.mass-scan.mass*error), xlinks.upper_bound(scan.mass+scan.mass*error) );
  
  map<double, ion> candidates;
  for(map<double, ion>::iterator it=submap.begin(); it!=submap.end(); ++it) {
    candidates.insert(pair<double, ion>(abs(scan.mass-(it->first)), it->second));
  }
  
  while(candidates.size() > top) {
    candidates.erase(--candidates.end());
  }
  
  xpsm x; x.scan = scan; x.candidates = candidates;
  return x;
}


vector<xpsm> match(map<double, ion> xlinks, map<double, ion> scans, double error = 0.00001, int top = 100) {

  vector<xpsm> xpsms;
  for(map<double, ion>::iterator it=scans.begin(); it!=scans.end(); ++it) {
  	xpsms.push_back(matchScan(xlinks, it->second, error, top));
  }
  
  return xpsms;
}


int main(const int argv, const char** args) {
  char sep = '\t';
  if(argv==4) { sep = args[3][0]; }
	
  map<double, ion> candidates = parse(args[1]);
//  for(map<double, ion>::iterator it=candidates.begin(); it!=candidates.end(); ++it) {
//    printf("%16.16f\n", it->first);
//  }

  cout << "============\n";

  map<double, ion> scans = parseMGF(args[2]);
  
  
  printf("scan_id%1$cscan_mass%1$cxlink_id%1$cxlink_mass%1$cerror\n", sep);
  
  vector<xpsm> xpsms = match(candidates, scans, 0.00001, 3);
  for(vector<xpsm>::iterator it=xpsms.begin(); it!=xpsms.end(); ++it) {
  	for(map<double, ion>::iterator x=it->candidates.begin(); x!=it->candidates.end(); ++x) {
  		//cout << it->scan.text << it->scan.mass << " " << x->first << " " << x->second.text << endl;
  		printf("%2$s%1$c%3$.8f%1$c%4$s%1$c%5$.8f%1$c%6$.8f\n", sep, it->scan.text.c_str(), it->scan.mass, x->second.text.c_str(), x->second.mass, x->first);
  	}
  }
  
  cout << "============\n";
//  xpsm matched = match(candidates, ion(13.23, "halo"), 10.1, 3);
//  for(map<double, ion>::iterator it=matched.candidates.begin(); it!=matched.candidates.end(); ++it) {
//    printf("%16.16f %s\n", it->first, it->second.text.c_str());
//  }
  return 0;
}
