/*
use AU.

5.03e-31 solar mass/kg

G = 38.83 AU^3/(SolarMass*yr^2)
*/

const int N = 2;//# of bodies
const double G = 39.478;//1.0;
const double tpi = 2*TMath::Pi();
const double t = 1.0/365.0;//pow(2*TMath::Pi(),-1);
const int T =  365*25; //days in 25 years
struct Particle{
	double mass;
	double R;
	double x,y,z;
	double xo,yo,zo;
	double xn,yn,zn;
	double vx,vy,vz;
	double vxo,vyo,vzo;
	double vxn,vyn,vzn;
	double ax,ay,az;
};

double Energy(double m, double vx, double vy, double vz){
	return .5*m*(vx*vx + vy*vy + vz*vz);
}

double AngMomentum(double m, double r1, double r2, double r3, double v1, double v2, double v3){
	double L,ihat,jhat,khat;
	
	ihat = r2*v3 - r3*v2;
	jhat = r3*v1 - r1*v3;
	khat = r1*v2 - r2*v1;
	
	L = m*pow(ihat*ihat + jhat*jhat + khat*khat,.5);
	
	return L;
}

double Distance(Particle *body1, Particle *body2){
	double xi,xj,yi,yj,zi,zj;
	double dx,dy,dz,R;
	xi = body1->x;
	yi = body1->y;
	zi = body1->z;
	
	xj = body2->x;
	yj = body2->y;
	zj = body2->z;
	
	dx = fabs(xi - xj);
	dy = fabs(yi - yj);
	dz = fabs(zi - zj);
	
	R = sqrt((dx*dx + dy*dy + dz*dz));
	
	return R;
}

void gravAccel(Particle *body1, Particle *body2)
{
	double dx,dy,dz,r,rp;
	double mi,mj,xi,xj,yi,yj,zi,zj;
	mi = body1->mass;
	xi = body1->x;
	yi = body1->y;
	zi = body1->z;
	
	mj = body2->mass;
	xj = body2->x;
	yj = body2->y;
	zj = body2->z;
	
	dx = xi - xj;
	dy = yi - yj;
	dz = zi - zj;
	//cout << "dx = " << dx << "dy = " << dy << "dz = " << dz << endl;
	
	r = pow((dx*dx + dy*dy + dz*dz),1.5);
	rp = sqrt((dx*dx + dy*dy + dz*dz));
	
	//cout << "r = " << rp << endl;
	dx = dx/r;
	dy = dy/r;
	dz = dz/r;
	
	body2->ax = body2->ax + G*mi*dx;
	body1->ax = body1->ax + G*mj*dx;
	
	body2->ay = body2->ay + G*mi*dy;
	body1->ay = body1->ay + G*mj*dy;
	
	body2->az = body2->az + G*mi*dz;
	body1->az = body1->az + G*mj*dz;
}

void nbody()
{
	 Particle *body = new Particle[N];
	 TGraphErrors* graphE = new TGraphErrors();
	 TGraphErrors* graphM = new TGraphErrors();
	 
	 TGraphErrors* graph1 = new TGraphErrors();
	 TGraphErrors* graph2 = new TGraphErrors();	 
	 TGraphErrors* graph3 = new TGraphErrors();
	 TGraphErrors* graph4 = new TGraphErrors();
	 TGraphErrors* graph5 = new TGraphErrors();
	 TGraphErrors* graph6 = new TGraphErrors();
	 TGraphErrors* graph7 = new TGraphErrors();
	 TGraphErrors* graph8 = new TGraphErrors();
	 TGraphErrors* graph9 = new TGraphErrors();
	 
	 
	 //Initial Values 
	 int a = 0;
	 int b = 0;
	 double time = 0;
	 double dEnergy,iEnergy,KinE,PotE;
	 double dMom, iMom;
	 string file = "2body.txt";
	 int length;
	 ifstream input;
	 cout <<" open file "<<file<<endl;
	 input.open(file.c_str());
	 cout << " file opend " << endl;
	 input >> length;
	 cout << " file size " << length << endl;
	
	 for(int i = 0;i < length; i++){		
		input >> body[b].mass;
		body[b].ax = 0;
		body[b].ay = 0;
		body[b].az = 0;
		//cout << body[b].mass << endl;	 
		
		input >> body[b].x;
		input >> body[b].y;
		input >> body[b].z;
		//cout << body[b].x << " " << body[b].y << " " << body[b].z << " " << endl;
		
	
		input >> body[b].vx;
		input >> body[b].vy;
		input >> body[b].vz;
		
		//convert au/day to au/yr
		body[b].vx = body[b].vx*365.0;///tpi;
		body[b].vy = body[b].vy*365.0;///tpi;
		body[b].vz = body[b].vz*365.0;///tpi;
		//cout << body[b].vx << " " << body[b].vy << " " << body[b].vz << " " << endl;
		b++;
	}

	graph1->SetPoint(0,body[1].x,body[1].y);
	graph2->SetPoint(0,body[2].x,body[2].y);
	graph3->SetPoint(0,body[3].x,body[3].y);
	graph4->SetPoint(0,body[4].x,body[4].y);
	graph5->SetPoint(0,body[5].x,body[5].y);
	graph6->SetPoint(0,body[6].x,body[6].y);
	graph7->SetPoint(0,body[7].x,body[7].y);
	graph8->SetPoint(0,body[8].x,body[8].y);
	graph9->SetPoint(0,body[9].x,body[9].y);
	
	//Initialize Energy & Angular Momentum
	iEnergy = 0.0;
	KinE = 0.0;
	PotE = 0.0;
	iMom = 0.0;
	for(int i = 0;i < N;i++){
		KinE = KinE + Energy(body[i].mass,body[i].vx,body[i].vy,body[i].vz);
		iMom = iMom + AngMomentum(body[i].mass,body[i].x,body[i].y,body[i].z,body[i].vx,body[i].vy,body[i].vz);
	}	
	
	for(int i = 0;i < N-1;i++){
		for(int j = i+1;j < N;j++){
			PotE = PotE + G*body[i].mass*body[j].mass/Distance(&body[i],&body[j]);
		}
	}	
	
	iEnergy = KinE - PotE;
	//cout << "E(0) = " << iEnergy << endl;
	//cout << "L(0) = " << iMom << endl;

	
	//Intitialize accerlations
	for(int i = 0;i < N-1;i++){
		for(int j = i+1;j < N;j++){
			gravAccel(&body[i],&body[j]);
		}
	}
	//cout<<"acc of body1: (" << body[0].ax << "," << body[0].ay << "," << body[0].az << ")" << endl;
	//cout<<"acc of body2: (" << body[1].ax << "," << body[1].ay << "," << body[1].az << ")" << endl;
	
	//1st step (1st deltat)
	for(int i = 0;i < N; i++){
		body[i].xo = body[i].x;
		body[i].yo = body[i].y;
		body[i].zo = body[i].z;
		
		body[i].vxo = body[i].vx;
		body[i].vyo = body[i].vy;
		body[i].vzo = body[i].vz;
		
		body[i].x = body[i].xo + t*body[i].vxo;
		body[i].vx = body[i].vxo + t*body[i].ax;
		
		body[i].y = body[i].yo + t*body[i].vyo;
		body[i].vy = body[i].vyo + t*body[i].ay;
		
		body[i].z = body[i].zo + t*body[i].vzo;
		body[i].vz = body[i].vzo + t*body[i].az;
		//cout << "v(1) = " << body[i].vx << " " << body[i].vy << " " << body[i].vz << " " << endl;
	}
	
	double Energy,Mom;
	//time after delt
	for(int k = 1;k <= T;k++){
		Energy = 0.0;
		Mom = 0.0;
		KinE = 0.0;
		PotE = 0.0;
	
		graph1->SetPoint(k,body[1].x,body[1].y);
		graph2->SetPoint(k,body[2].x,body[2].y);
		graph3->SetPoint(k,body[3].x,body[3].y);
		graph4->SetPoint(k,body[4].x,body[4].y);
		graph5->SetPoint(k,body[5].x,body[5].y);
		graph6->SetPoint(k,body[6].x,body[6].y);
		graph7->SetPoint(k,body[7].x,body[7].y);
		graph8->SetPoint(k,body[8].x,body[8].y);
		
		graph9->SetPoint(k,body[9].x,body[9].y);

		//Find Energy & Momentum
		for(int i = 0;i < N;i++){
			//for(int j = i+1;j < N;j++){
				//PotE = PotE + G*body[i].mass*body[j].mass/Distance(&body[i],&body[j]);
			//}
			KinE = KinE + Energy(body[i].mass,body[i].vx,body[i].vy,body[i].vz);
			Mom = Mom + AngMomentum(body[i].mass,body[i].x,body[i].y,body[i].z,body[i].vx,body[i].vy,body[i].vz);
		}	
		
		for(int i = 0;i < N-1;i++){
			for(int j = i+1;j < N;j++){
				PotE = PotE + G*body[i].mass*body[j].mass/Distance(&body[i],&body[j]);
			}
		}	
		Energy = KinE - PotE;
		//Find dE
		//cout << "E(" << k << ") = " << Energy << endl;
		dEnergy = (iEnergy-Energy)/iEnergy;
		//cout << "dE = " << dEnergy << endl;
		
		//Find dL
		//cout << "L(" << k << ") = " << Mom << endl;
		dMom = (Mom - iMom)/iMom;
		//cout << "dL = " << dMom << endl;
		
		//reintialize acceleration for the time.
		for(int i = 0;i < N;i++){
			body[i].ax = 0;
			body[i].ay = 0;
			body[i].az = 0;
		}
		//cout << "v("<< k << ") = (" << body[1].vx << "," << body[1].vy << "," << body[1].vz << ")" << endl;
		
		//run acceleration routine
		for(int i = 0;i < N-1;i++){
			for(int j = i+1;j < N;j++){
				gravAccel(&body[i],&body[j]);
			}
		}
	//cout<<"acc of body1: (" << body[0].ax << "," << body[0].ay << "," << body[0].az << ")" << endl;
	//out<<"acc of body2: (" << body[1].ax << "," << body[1].ay << "," << body[1].az << ")" << endl;
	
		//leap frog
		for(int i = 0;i < N; i++){
			body[i].xn = body[i].xo + 2.0*t*body[i].vx;
			body[i].vxn = body[i].vxo + 2.0*t*body[i].ax;
			
			body[i].yn = body[i].yo + 2.0*t*body[i].vy;
			body[i].vyn = body[i].vyo + 2.0*t*body[i].ay;
			
			body[i].zn = body[i].zo + 2.0*t*body[i].vz;
			body[i].vzn = body[i].vzo + 2.0*t*body[i].az;
			//cout << body[i].mass << endl;
			//cout << body[i].xn << " " << body[i].yn << " " << body[i].zn << " " << endl;
			//cout << body[i].vxn << " " << body[i].vyn << " " << body[i].vzn << " " << endl;
		}
		time = time+t*(365.0);
		
		//Plot Points for Energy and Mom
		graphE->SetPoint(k-1,time/365.0,dEnergy);
		graphM->SetPoint(k-1,time/365.0,dMom);
		
		
		//reinitialize
		for(int i=0;i<N;i++){
			body[i].xo = body[i].x;
			body[i].yo = body[i].y;
			body[i].zo = body[i].z;
			
			body[i].vxo = body[i].vx;
			body[i].vyo = body[i].vy;
			body[i].vzo = body[i].vz;
			
			body[i].x = body[i].xn;
			body[i].y = body[i].yn;
			body[i].z = body[i].zn;
			
			body[i].vx = body[i].vxn;
			body[i].vy = body[i].vyn;
			body[i].vz = body[i].vzn;
			
		}
	}
	
	cout << "time: " << time << endl;
	//Plot Energy
	TCanvas *c1 = new TCanvas("c1","Assignment #2",200,10,700,500);
   //c1->SetLogy();
	//c1->SetLogx();

   graphE->Draw("AP");
   graphE->SetMarkerColor(kRed);
   graphE->SetMarkerStyle(3);
   graphE->SetMarkerSize(.2);
   graphE->SetLineColor(kRed);
   graphE->SetLineWidth(0.1);
   
   graphE->GetXaxis()->SetTitle("Years");
	graphE->GetYaxis()->SetTitle("#DeltaE/E_{0}");
	graphE->GetXaxis()->CenterTitle();
	graphE->GetYaxis()->CenterTitle();
	//graphE->GetYaxis()->SetRangeUser(1e-7,1e-5);
  
   //Plot Mom
	TCanvas *c2 = new TCanvas("c2","Assignment #2",200,10,700,500);
   //c2->SetLogy();
	//c1->SetLogx();

   graphM->Draw("AP");
   graphM->SetMarkerColor(kRed);
   graphM->SetMarkerStyle(3);
   graphM->SetMarkerSize(.2);
   graphM->SetLineColor(kRed);
   graphM->SetLineWidth(0.1);
   
   
   graphM->GetXaxis()->SetTitle("Years");
   graphM->GetYaxis()->SetTitle("#DeltaL/L_{0}");
	graphM->GetXaxis()->CenterTitle();
	graphM->GetYaxis()->CenterTitle();
	//graphM->GetYaxis()->SetRangeUser(1e-7,1e-5);
  
 /* 
  	TMultiGraph *mgr = new TMultiGraph();
  	mgr->Add(graph1);
   mgr->Add(graph2);
   mgr->Add(graph3);
   mgr->Add(graph4);
   mgr->Add(graph5);
   mgr->Add(graph6);
   //mgr->Add(graph7);
   //mgr->Add(graph8);
   //mgr->Add(graph9);
   mgr->Draw("AP");
 	mgr->GetXaxis()->SetTitle("AU");
	mgr->GetYaxis()->SetTitle("AU");
	mgr->GetXaxis()->CenterTitle();
	mgr->GetYaxis()->CenterTitle();
 	graph5->SetMarkerColor(kRed);
*/ 
   delete [] body;
   
}

