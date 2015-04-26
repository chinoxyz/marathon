#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <set>
#include <queue>
#include <stack>
#include <string>
#include <sstream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cassert>

using namespace std;

#define ri(X) scanf("%d", &(X))
#define pi(X) printf("%d", (X))
#define mp(X,Y) make_pair(X,Y)
#define pb(X) push_back(X)
#define lint long long
#define pii pair<int,int>
#define inf 1e9
#define linf 1e18
#define X first
#define Y second
#define all(X) (X).begin(),(X).end()
#define uni(X) X.erase(unique(X.begin(), X.end()), X.end());
double INF = 1e100;
double EPS = 1e-12;






string itos(int x){
  stringstream ss;
  ss << x;
  return ss.str();
}


const double pi = acos(-1);
class PT{
  public:
  double x, y;
  int i;
  PT(int xx = 0, int yy = 0, int ii = 0){
    x = xx;
    y = yy;
    i = ii;
  }
  bool operator <(const PT &p)  const { return x < p.x || (x == p.x && y < p.y); }
  bool operator ==(const PT &p)  const { return x == p.x && y == p.y; }
  PT operator + (const PT &p)  const { return PT(x+p.x, y+p.y); }
  PT operator - (const PT &p)  const { return PT(x-p.x, y-p.y); }
  PT operator * (double c)     const { return PT(x*c,   y*c  ); }
  PT operator / (double c)     const { return PT(x/c,   y/c  ); }
};
double dot(PT p, PT q)     { return p.x*q.x+p.y*q.y; }
double dist2(PT p, PT q)   { return dot(p-q,p-q); }
double dist(PT p, PT q)   { return sqrt(dist2(p,q)); }
double cross(PT p, PT q)   { return p.x*q.y-p.y*q.x; }

double angle0(PT p){
  return atan2(p.x,p.y);
}

ostream &operator<<(ostream &os, const PT &p) {
  os << "(" << p.x << "," << p.y << ")"; 
}

// O(ni)
double ComputeArea(const vector<PT> &p) {
  double area = 0;
  for(int i = 0; i < p.size(); i++) {
    int j = (i+1) % p.size();
    area += p[i].x*p[j].y - p[j].x*p[i].y;
  }
  return fabs(area / 2.0);
}
// O(1)
bool LinesParallel(PT a, PT b, PT c, PT d) { 
  return fabs(cross(b-a, c-d)) < EPS; 
}

// O(1)
bool LinesCollinear(PT a, PT b, PT c, PT d) { 
  return LinesParallel(a, b, c, d)
      && fabs(cross(a-b, a-c)) < EPS
      && fabs(cross(c-d, c-a)) < EPS; 
}

// O(1)
bool SegmentsIntersect(PT a, PT b, PT c, PT d) {
  if (LinesCollinear(a, b, c, d)) {
    if (dist2(a, c) < EPS || dist2(a, d) < EPS ||
      dist2(b, c) < EPS || dist2(b, d) < EPS) return true;
    if (dot(c-a, c-b) > 0 && dot(d-a, d-b) > 0 && dot(c-b, d-b) > 0)
      return false;
    return true;
  }
  if (cross(d-a, b-a) * cross(c-a, b-a) > 0) return false;
  if (cross(a-c, d-c) * cross(b-c, d-c) > 0) return false;
  return true;
}

// O(ni)
bool hasArea(const vector<PT> &p){
  bool allcol = 1;
  for(int i = 0; i < p.size()-1 && allcol; i++){
    if(!LinesParallel(p[i], p[i+1], p[i+1], p[(i+2)%p.size()])){
      allcol = 0;
    }
  }
  if(allcol) return 0;
  return 1;
}
// O(ni*ni)
bool IsSimple(const vector<PT> &p) {
  
  
  for (int i = 0; i < p.size(); i++) {
    for (int k = i+1; k < p.size(); k++) {
      int j = (i+1) % p.size();
      int l = (k+1) % p.size();
      if (i == l || j == k) continue;
      if (SegmentsIntersect(p[i], p[j], p[k], p[l])){
      
        for(int w = 0; w < p.size(); w++){
          cerr<< p[w] << " ";
        }cerr << endl;
      
        cerr << p[i] << p[j] << p[k] << p[l] << endl;
        return false;
      }
    }
  }
  return true;
}
vector<PT> ve;
vector<string> res;
int np;
int n;


class poly{
public:
  vector<int> vi;
  double cost;
  poly(vector<int> vv){
    cost = -1;
    vi = vv;
  }
  
  // O(C(Ri))
  double getCost(){
    if(cost >= 0){
      return cost;
    }
    vector<PT> vq;
    for(int i = 0; i < vi.size(); i++){
      vq.pb(ve[vi[i]]);
    }
    if(!hasArea(vq)){
      cost = inf;
      return cost;
    }
    cost = ComputeArea(vq);
    return cost;
  }
  // O(npi)
  string tostring(){
    string rp;
    for(int i = 0; i < vi.size(); i++){
      rp+=itos(vi[i])+ " ";
    }
    return rp;
  }
  
};


// O(ni lg ni)
poly getpol1(vector<int> vi, int s=0, int sk=0){
  assert(vi.size() >= 2);
  if(vi.size() == 3) return vi;
  //assert(np == ve.size());
  assert(np == vi.size());
  PT prom = PT(0,0), vprom;
  for(int i = 0; i < vi.size(); i++){
    prom=prom+ve[vi[i]];
  }
  prom = prom/vi.size();
  cout << prom << endl;
  vprom = prom/100000;
  
  PT ic = ve[vi[s]];//PT(pi/12,pi/15);
  cout << "CENTER: "<<ic << endl;
  
  vector<pair<pair<double,double>, int> > vs;
  for(int i = 0; i < vi.size(); i++){
    if(vi[i] == vi[s]){
      continue;
    }
    double ang = angle0(ve[vi[i]]-ic);
    double di = dist2(ic, ve[vi[i]]);
    vs.pb(mp(mp(ang,di),vi[i]));
  }
  sort(all(vs));
  vi[0] = vi[s];
  
  int ki = 0;
  for(int i = 0; i < vs.size()-1; i++){
    
    if(fabs(vs[i].X.X-vs[i+1].X.X) >= pi){
      cerr << vs[i].X.X << " " << vs[i+1].X.X << endl;
      ki = i+1;
      cerr << "ki " << ki << endl;
    }
  }
  
  
  for(int i = 0; i < vs.size(); i++){
    //vi[i] = vs[(i)%vs.size()].Y;
    vi[i+1] = vs[(i+ki)%vs.size()].Y;
  }
  
  for(int i = 0; i < vi.size(); i++){
    cout << ve[vi[i]] << " ";
  }cout << endl;
  
  vector<pair<double, int> > vc;
  int si = vi.size()-1;
  for(int k = vi.size()-1; k >= 0; k--){
    if(LinesParallel(ve[vi[0]], ve[vi[vi.size()-1]], ve[vi[(k+1)%vi.size()]], ve[vi[k]])){
      cout << ve[vi[0]] << ve[vi[vi.size()-1]]<< ve[vi[(k+1)%vi.size()]] << ve[vi[k]]<< " ";
      double di = dist2(ve[vi[0]], ve[vi[k]]);
      vc.pb(mp(-di,vi[k]));
      si = k;
    }else{
      break;
    }
  }cout << endl;
  

    
  
  sort(all(vc));
  
  cout << "si:"<<si << endl;
  for(int i = 0; i < vc.size(); i++){
    cout << ve[vc[i].Y] << " ";
  }cout << endl;
  for(int i = 0; i < vc.size(); i++){
    vi[si+i] = vc[i].Y;
  }
  
  
  
  
  
  
  return poly(vi);
}








class SmallPolygons{
  
  
  public:
  void choosePolygons(vector <int> points, int N){
    ve.clear();
    np = points.size();
    n = N;
    for(int i = 0; i < np; i+=2){
      ve.pb(PT(points[i],points[i+1],i/2));
    }
    np/=2;
    
    
    vector<int> vi;
    for(int i = 0; i < np; i++){
      vi.pb(i);
    }
    for(int i = 0; i < np; i++){
      poly w = getpol1(vi,i);
      
      vector<PT> vw;
      for(int k = 0; k < np; k++){
        vw.pb(ve[w.vi[k]]);
      }
      
      if(!IsSimple(vw)){
        cout << "hola: " << i << endl;
        exit(0);
      }
      cout << np << " " << ve.size() << "cic" << i << endl;
    }
    
  }

};


int main(){
  
  
  for(int k = 0; k < 1000000; k++){
    int n = 6;
    int mn = 100;
    vector<pii> ve;
    for(int i = 0; i < n; i++){
      int x,y;
      
      x = rand()%mn;
      y = rand()%mn;
      
      ve.pb(mp(x,y));
    }
    sort(all(ve));
    uni(ve);
    vector<int> vi;
    for(int i = 0; i < ve.size(); i++){
      vi.pb(ve[i].X);
      vi.pb(ve[i].Y);
    }
    for(int i = 0; i < ve.size(); i++){
      printf("[%d %d] ", ve[i].X, ve[i].Y);
    }cout << endl;
    
    SmallPolygons().choosePolygons(vi, 1);
  }
  return 0;
}

