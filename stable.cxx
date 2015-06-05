/*
use AU.

5.03e-31 solar mass/kg

G = 39.478 AU^3/(SolarMass*yr^2)


*/

const int N = 3;//# of bodies
const double G = 1.0;
const double tpi = 2.0*TMath::Pi();
const double t = tpi/365.0;//pow(2.0*TMath::Pi(),-1);
const int T = 200;

struct Particle{
	double mass;
	double x,y,z;
	double xo,yo,zo;
	double xn,yn,zn;
	double vx,vy,vz;
	double vxo,vyo,vzo;
	double vxn,vyn,vzn;
	double ax,ay,az;
};

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
	body1->ax = body1->ax - G*mj*dx;
	
	body2->ay = body2->ay + G*mi*dy;
	body1->ay = body1->ay - G*mj*dy;
	
	body2->az = body2->az + G*mi*dz;
	body1->az = body1->az - G*mj*dz;
}

void stable()
{
	 Particle *body = new Particle[N];
	 TGraphErrors* graph1 = new TGraphErrors();
	 TGraphErrors* graph2 = new TGraphErrors();
	 TGraphErrors* graph3 = new TGraphErrors();
	 
	 //Initial Values 
	 int a = 0;
	 int b = 0;
	 int time = 0;
	 
	 string file = "3body.txt";
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
		
		input >> body[b].vx;
		input >> body[b].vy;
		input >> body[b].vz;
		
		//body[b].vx = body[b].vx/tpi;
		//body[b].vy = body[b].vy/tpi;
		//body[b].vz = body[b].vz/tpi;
		
		//cout << body[b].x << " " << body[b].y << " " << body[b].z << " " << endl;
		//cout << body[b].vx << " " << body[b].vy << " " << body[b].vz << " " << endl;
		b++;
	 }

		graph1->SetPoint(0,body[0].x,body[0].y);
		graph2->SetPoint(0,body[1].x,body[1].y);
		graph3->SetPoint(0,body[2].x,body[2].y);
		
		//Intitialize accerlations
		for(int i = 0;i < N;i++){
			for(int j = i+1;j < N;j++){
				gravAccel(&body[i],&body[j]);
			}
		}
		//cout<<"acc of body1: (" << body[0].ax << "," << body[0].ay << "," << body[0].az << ")" << endl;
		//cout<<"acc of body2: (" << body[1].ax << "," << body[1].ay << "," << body[1].az << ")" << endl;
		
		//1st step (1st year)
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
		
//==================================================================================================================		
		//time after 1st year
		for(int k = 1;k <= T;k++){
			
			graph1->SetPoint(k,body[0].x,body[0].y);
			graph2->SetPoint(k,body[1].x,body[1].y);
			graph3->SetPoint(k,body[2].x,body[2].y);
			
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
		//cout<<"acc of body2: (" << body[1].ax << "," << body[1].ay << "," << body[1].az << ")" << endl;
		
		//leap frog
		for(int i = 0;i < N; i++){
			body[i].xn = body[i].xo + 2.0*t*body[i].vx;
			body[i].vxn = body[i].vxo + 2.0*t*body[i].ax;
			
			body[i].yn = body[i].yo + 2.0*t*body[i].vy;
			body[i].vyn = body[i].vyo + 2.0*t*body[i].ay;
			
			body[i].zn = body[i].zo + 2.0*t*body[i].vz;
			body[i].vzn = body[i].vzo + 2.0*t*body[i].az;
			//cout << body[i].mass << endl;
			//cout << "position: " <<  body[i].xn << " " << body[i].yn << " " << body[i].zn << " " << endl;
			//cout << "velocity: " << body[i].vxn << " " << body[i].vyn << " " << body[i].vzn << " " << endl;
		}
		
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
//==============================================================================================================	
	TCanvas *c1 = new TCanvas("c1","Assignment #2",200,10,700,500);
   //c1->SetLogy();
   
   //graph1->Draw("AP");
   graph1->SetMarkerColor(kRed);
   graph1->SetMarkerStyle(3);
   graph1->SetMarkerSize(1);
   graph1->SetLineColor(kRed);
   graph1->SetLineWidth(0.1);
   
   graph2->SetMarkerColor(kBlue);
   graph2->SetMarkerStyle(3);
   graph2->SetMarkerSize(1);
   graph2->SetLineColor(kBlue);
   graph2->SetLineWidth(0.1);
   
   graph3->SetMarkerColor(kViolet);
   graph3->SetMarkerStyle(3);
   graph3->SetMarkerSize(1);
   graph3->SetLineColor(kViolet);
   graph3->SetLineWidth(0.1);
   
   TMultiGraph *mgr = new TMultiGraph();
   mgr->Add(graph1);
   mgr->Add(graph2);
   mgr->Add(graph3);
   //mgr->SetTitle("Time = 365 days");
   mgr->Draw("AP");
   mgr->GetXaxis()->SetTitle("AU");
	mgr->GetYaxis()->SetTitle("AU");
	mgr->GetXaxis()->CenterTitle();
	mgr->GetYaxis()->CenterTitle();
   delete [] body;
   
}