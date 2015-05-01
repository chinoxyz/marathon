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


vector<PT> hull(vector<PT> P, bool b = false){
  assert(P.size() >= 3);
	int n = P.size(), k = 0; 
	vector<PT> H(2*n); 
	// Sort PTs lexicographically 
	sort(P.begin(), P.end()); 
	// Build lower hull 
	for (int i = 0; i < n; i++) {
	  if(b == 0){
		  while (k >= 2 && cross(H[k-1]-H[k-2], P[i]-H[k-2]) < 0) k--; 
		}else{
		  while (k >= 2 && cross(H[k-1]-H[k-2], P[i]-H[k-2]) > 0) k--; 
		}
		H[k++] = P[i]; 
	}
	H.resize(k);
	return H; 
}

vector<PT> convex_hull(vector<PT> P){
  assert(P.size() >= 3);
	int n = P.size(), k = 0; 
	vector<PT> H(2*n); 
	// Sort PTs lexicographically 
	sort(P.begin(), P.end()); 
	// Build lower hull 
	for (int i = 0; i < n; i++) { 
		while (k >= 2 && cross(H[k-1]-H[k-2], P[i]-H[k-2]) < 0) k--; 
		H[k++] = P[i]; 
	} 
	// Build upper hull 
	for (int i = n-2, t = k+1; i >= 0; i--) { 
		while (k >= t && cross(H[k-1]-H[k-2], P[i]-H[k-2]) < 0) k--; 
		H[k++] = P[i]; 
	} 
	H.resize(k);
	H.pop_back();
	return H; 
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
  poly(){
    cost = inf;
  }
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
    if(!hasArea(vq) || !IsSimple(vq)){
      //cerr << "SHiT" << endl;
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

poly getpol1(vector<int>, int, int);
poly getpol(vector<int> vi);

vector<int> reorder(PT center, vector<int> vi){
  int ki = 0;
  for(int i = 0; i < vi.size(); i++){
    if(ve[vi[i]] == center){
      ki = i;
    }
  }
  cout << center <<" " <<  ki << endl;
  if(ki == 0) return vi;
  rotate(vi.begin(), vi.begin()+ki, vi.end());
  
  assert(ve[vi[0]] == center);
  
  return vi;
}


class solution{
  public:
  vector<poly> vp;
  double cost;
  solution(){
    cost = inf;
  }
  solution(vector<poly> vv){
    vp = vv;
    cost = -1;
  }
  
  // O(C(R)) 
  double getCost(){
    if(cost >= 0){
      return cost;
    }
    cost = 0;
    for(int i = 0; i < vp.size(); i++){
      cost+=vp[i].getCost();
    }
    return cost;
  }
  
  // O(np)
  vector<string> tovs(){
    vector<string> vs;
    for(int i = 0; i < vp.size(); i++){
      vs.pb(vp[i].tostring());
    }
    return vs;
  }
  
  
  
  
  
  
  // O(ni lg ni + C(R));
  void split(int k, int w){
    ;if(vp.size() >= n) return;
    assert(0 <= k && k < vp.size());
    if(vp[k].vi.size() <= 6) return;
    //cerr <<"Hola:  "<< k << " " << vp.size() <<" " << endl;
    int rot = rand()%2;
    vector<int>&vi = vp[k].vi;
    vector<pair<pair<int,int>,int> > vs;
    for(int i = 0; i < vi.size(); i++){
      if(rot == 0){
        vs.pb(mp(mp(ve[vi[i]].x,ve[vi[i]].y),vi[i]));
      }else{
        vs.pb(mp(mp(ve[vi[i]].y,ve[vi[i]].x),vi[i]));
      }
    }
    sort(all(vs));
    vector<int> rp[2];
    for(int i = 0; i < vs.size(); i++){
      rp[i< w].pb(vs[i].Y);//
    }
    vp[k] = getpol(rp[0]);
    vp.pb(getpol(rp[1]));
    cost = -1;
  }
  
  // O(ni lg ni + c(R))
  void splitr(int k){
    assert(0 <= k && k < vp.size());
    int siz = vp[k].vi.size();
    if(siz <= 8) return;
    split(k, 3+(rand()%(siz-6)));
  }
  
  // O(ni lg ni + C(R))
  void star(int k){
    assert(0 <= k && k < vp.size());
    vp[k] = getpol(vp[k].vi);
    cost = -1;
  }
};

// O(ni lg ni)
poly getpol1(vector<int> vi, int s=0){
  assert(vi.size() >= 2);
  if(vi.size() == 3) return vi;
  PT ic = ve[vi[s]];
  cout << "center:"<<ic << endl;
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
      ki = i+1;
    }
  }
  for(int i = 0; i < vs.size(); i++){
    vi[i+1] = vs[(i+ki)%vs.size()].Y;
  }
  vector<pair<double, int> > vc;
  int si = vi.size()-1;
  for(int k = vi.size()-1; k >= 0; k--){
    if(LinesParallel(ve[vi[0]], ve[vi[vi.size()-1]], ve[vi[(k+1)%vi.size()]], ve[vi[k]])){
      double di = dist2(ve[vi[0]], ve[vi[k]]);
      vc.pb(mp(-di,vi[k]));
      si = k;
    }else{
      break;
    }
  }
  
  sort(all(vc));
  for(int i = 0; i < vc.size(); i++){
    vi[si+i] = vc[i].Y;
  }
  
  return poly(reorder(ic,vi));
}

// O(ni lg ni)
poly getpol(vector<int> vi){
  return getpol1(vi, rand()%vi.size());
}


// O(ni lg ni + C(R));
poly getpol2(vector<int> vi){
  if(vi.size() <= 6) return getpol(vi);
  int w = vi.size()/2;
  //cerr <<"Hola:  "<< k << " " << vp.size() <<" " << endl;
  int rot = 1;//rand()%2;
  vector<pair<pair<int,int>,int> > vs;
  for(int i = 0; i < vi.size(); i++){
    if(rot == 0){
      vs.pb(mp(mp(ve[vi[i]].x,ve[vi[i]].y),vi[i]));
    }else{
      vs.pb(mp(mp(ve[vi[i]].y,ve[vi[i]].x),vi[i]));
    }
  }
  sort(all(vs));
  vector<int> rp[2];
  vector<PT> vh[2];
  for(int i = 0; i < vs.size(); i++){
    rp[i< w?0:1].pb(vs[i].Y);//
    vh[i< w?0:1].pb(ve[vs[i].Y]);//
  }
  //cout << w << endl;
  vh[0] = hull(vh[0],1);
  vh[1] = hull(vh[1],0);
  
  for(int k = 0; k < 2; k++){
    for(int i = 0; i < rp[k].size(); i++){
      cout << ve[rp[k][i]] << " ";
    }cout << endl;
    for(int i = 0; i < vh[k].size(); i++){
      cout << vh[k][i] << " ";
    }cout << endl;
    cout << endl;
  }
  
  
  rp[0] = reorder(vh[0][1], rp[0]);
  rp[1] = reorder(vh[1][1], rp[1]);
  
  
  for(int k = 0; k < 2; k++){
    for(int i = 0; i < rp[k].size(); i++){
      cout << ve[rp[k][i]] << " ";
    }cout << endl;
    for(int i = 0; i < vh[k].size(); i++){
      cout << vh[k][i] << " ";
    }cout << endl;
    cout << endl;
  }
  
  
  rp[0] = getpol1(rp[0],0).vi;
  rp[1] = getpol1(rp[1],0).vi;
  
  reverse(all(rp[1]));
  
  rp[0] = reorder(vh[0][0], rp[0]);
  rp[1] = reorder(vh[1][0], rp[1]);
  
  
  
  
  
  
  
  for(int k = 0; k < 2; k++){
    for(int i = 0; i < rp[k].size(); i++){
      cout << ve[rp[k][i]] << " ";
    }cout << endl;
    for(int i = 0; i < vh[k].size(); i++){
      cout << vh[k][i] << " ";
    }cout << endl;
    cout << endl;
  }
  
  
  int k = 0;
  
  for(int i = 0; i < rp[0].size(); i++){
    vi[k++] = rp[0][i];
  }
  for(int i = 1; i < rp[1].size(); i++){
    vi[k++] = rp[1][i];
  }
  vi[k++] = rp[1][0];
  
  
  if(!SegmentsIntersect(ve[rp[0][0]],ve[rp[1][0]], ve[rp[0][1]],ve[rp[1][1]]) ||
     SegmentsIntersect(ve[rp[0][0]],ve[rp[1][0]], ve[rp[0][0]]+PT(-1,0),ve[rp[0][1]]) ||
     SegmentsIntersect(ve[rp[0][1]],ve[rp[1][1]], ve[rp[0][0]],ve[rp[0][1]]+PT(-1,0))  
     ){
    cout << "HOLA!" << endl;
    return getpol1(vi);
  }
  cout << "YYYYYYYYYYYYYYYYYYEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEY" << endl;
  return vi;
  exit(0);
  
  /*
  for(int i = 0; i < pv[0].vi.size(); i++){
    cout  << ve[pv[0].vi[i]] << " ";
  }cout << endl;
  for(int i = 0; i < pv[1].vi.size(); i++){
    cout  << ve[pv[1].vi[i]] << " ";
  }cout << endl;
  int k = 0;
  
  for(int i = 0; i < pv[0].vi.size(); i++){
    vi[k++] = pv[0].vi[i];
  }
  for(int i = 1; i < pv[1].vi.size(); i++){
    vi[k++] = pv[1].vi[i];
  }
  vi[k++] = pv[1].vi[0];
  
  
  
  
  /*if(!SegmentsIntersect(pv[0].vi[0],pv[1].vi[0], pv[0].vi[1],pv[1].vi[1])){
    cout << "HEEY" << endl;
    return getpol1(vi);
  }
  vi[k++] = pv[0].vi[0];
  for(int i = 1; i < pv[1].vi.size(); i++){
    vi[k++] = pv[1].vi[i];
  }
  for(int i = 1; i < pv[0].vi.size(); i++){
    vi[k++] = pv[0].vi[i];
  }
  vi[k++] = pv[1].vi[0];
  */
  return getpol(vi);
}
  


// O(ni lg ni + C(R));
poly getpol3(vector<int> vi){
  if(vi.size() <= 6) return getpol(vi);
  set<int> se;
  vector<PT> vw,vc1,vc2,vr;
  vector<int> vi2;
  for(int i = 0; i < vi.size(); i++){
    vw.pb(ve[vi[i]]);
    se.insert(vi[i]);
  }
  vc1 = convex_hull(vw);
  
  for(int i = 0; i < vc1.size(); i++){
    se.erase(vc1[i].i);
  }
  if(se.size() <= 4){
    return getpol(vi);
  }
  
  for(set<int>::iterator it = se.begin(); it != se.end(); it++){
    vr.pb(ve[*it]);
  }
  vc2 = convex_hull(vr);
  
  vector<int> res;
  
  
  
  
  
  
  for(int i = 0; i < vc2.size(); i++){
    vi2.pb(vc2[i].i);
    se.erase(vc2[i].i);
  }
  for(set<int>::iterator it = se.begin(); it != se.end(); it++){
    vi2.pb(ve[*it].i);
  }
  vi2 = getpol1(vi2,0).vi;
  
  cout  << "CONVEXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX:"<<endl;
  for(int i = 0; i < vc1.size(); i++){
    cout << vc1[i] << " ";
  }cout << endl;
  for(int i = 0; i < vi2.size(); i++){
    cout << ve[vi2[i]] << " ";
  }cout << endl;
  
  //cout << LinesParallel(ve[vi2[0]], ve[vi2[1]], ve[vi2[1]], ve[vi2[2]]) << endl;
  
  if(!SegmentsIntersect(ve[vi2[0]], vc1.back(), ve[vi2[1]], vc1[0]) || LinesParallel(ve[vi2[0]], ve[vi2[1]], ve[vi2[1]], ve[vi2[2]]) ){
    return getpol(vi);
  }
  
  for(int i = 1; i < vi2.size(); i++){
    res.pb(vi2[i]);
  }
  res.pb(vi2[0]);
  for(int i = 0; i < vc1.size(); i++){
    res.pb(vc1[i].i);
  }
  assert(res.size() == vi.size());
  return res;
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
    
    poly w = getpol2(vi);
      
    
    
    vector<PT> vw;
    for(int k = 0; k < np; k++){
      vw.pb(ve[w.vi[k]]);
    }
    
    if(!IsSimple(vw)){
      cout << "hola: no es simple"  << endl;
      exit(0);
    }
    cout << np << " " << ve.size() << endl;
    /*
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
    }*/
    
  }

};


int main(){
  
  
  for(int k = 0; k < 1000000; k++){
    int n = 10;
    int mn = 30;
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

