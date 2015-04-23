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
  int x, y, i;
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
double angle(PT o, PT p, PT q){
  return ((cross(p-o,q-o)<0?-1:1)*acos(dot(p-o,q-o)/(dist(o,p)*dist(o,q))));
}
double angle2(PT o, PT p, PT q){
  return asin(cross(p-o,q-o)/(dist(o,p)*dist(o,q)));
}
ostream &operator<<(ostream &os, const PT &p) {
  os << "(" << p.x << "," << p.y << ")"; 
}
double ComputeArea(const vector<PT> &p) {
  double area = 0;
  for(int i = 0; i < p.size(); i++) {
    int j = (i+1) % p.size();
    area += p[i].x*p[j].y - p[j].x*p[i].y;
  }
  return fabs(area / 2.0);
}
bool LinesParallel(PT a, PT b, PT c, PT d) { 
  return fabs(cross(b-a, c-d)) < EPS; 
}

bool LinesCollinear(PT a, PT b, PT c, PT d) { 
  return LinesParallel(a, b, c, d)
      && fabs(cross(a-b, a-c)) < EPS
      && fabs(cross(c-d, c-a)) < EPS; 
}
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
bool IsSimple(const vector<PT> &p) {
  for (int i = 0; i < p.size(); i++) {
    for (int k = i+1; k < p.size(); k++) {
      int j = (i+1) % p.size();
      int l = (k+1) % p.size();
      if (i == l || j == k) continue;
      if (SegmentsIntersect(p[i], p[j], p[k], p[l])){
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
  double getCost(){
    if(cost >= 0){
      return cost;
    }
    vector<PT> vq;
    for(int i = 0; i < vi.size(); i++){
      vq.pb(ve[vi[i]]);
    }
    if(!IsSimple(vq)){
      cost = inf;
      return cost;
    }
    cost = ComputeArea(vq);
    return cost;
  }
  string tostring(){
    string rp;
    for(int i = 0; i < vi.size(); i++){
      rp+=itos(vi[i])+ " ";
    }
    return rp;
  }
  
};


class solution{
  public:
  vector<poly> vp;
  double cost;
  solution(){
    cost = -1;
  }
  solution(vector<poly> vv){
    vp = vv;
    cost = -1;
  }
  double getCost(){
    if(cost >= 0){
      return cost;
    }
    cost = 0;
    for(int i = 0; i < vp.size(); i++){
      vp[i].getCost();
    }
    return cost;
  }
  vector<string> tovs(){
    vector<string> vs;
    for(int i = 0; i < vp.size(); i++){
      vs.pb(vp[i].tostring());
    }
    return vs;
  }

  
};



poly getpol1(vector<int> vi, int s=0,PT pa = PT(-1000,0)){
  assert(vi.size() >= 3);
  if(vi.size() == 3) return vi;
  vector<PT> vp;
  for(int i = 0;i < vi.size(); i++){
    vp.pb(ve[vi[i]]);
  }
  sort(all(vp));
  PT ic = vp[s];
  PT ir = ic+pa;
  vector<pair<pair<double,double>, int> > vs;
  for(int i = 0; i < vi.size(); i++){
    
    if(ve[vi[i]] == ic){
      vs.pb(mp(mp(-1,-1),vi[i]));
      continue;
    } 
    double ang = angle(ic, ve[vi[i]], ir);
    double di = dist2(ic, ve[vi[i]]);
    vs.pb(mp(mp(ang,di),vi[i]));
  }
  sort(all(vs));
  for(int i = 0; i < vs.size(); i++){
    vi[i] = vs[i].Y;
  }
  return poly(vi);
}

poly getpol(vector<int> vi){
  return getpol1(vi, rand()%vi.size(), PT((rand()%1000)-500, (rand()%1000)-500));
}



solution getsol1(){
  vector<poly> res;
  vector<int> vi;
  for(int i = 0; i < np; i++){
    vi.pb(i);
  }
  res.pb(getpol1(vi));
  return solution(res);
}
/*
vector<vector<int> > getsol2(){
  //return getsol1();
  if(n == 1) return getsol1();
  
  vector<vector<int> > res;
  vector<PT> vp = ve;
  sort(all(vp));
  
  vector<int> vb;
  for(int i = 0; i < n-1; i++){
    vb.pb(rand()%np);
  }vb.pb(np-1);
  
  sort(all(vb));
  uni(vb);
  reverse(all(vb));
  
  
  for(int i = 0; i < 2; i++){
    res.pb(vector<int>());
  }
  
  for(int i = 0; i < np; i++){
    res[(i<np/2?0:1)].pb(vp[i].i);//
    /*vi.pb(i);//vp[i].i
    fprintf(stderr,"%d ",vp[i].i);
    if(!vb.empty() && vb.back() == i && vi.size() >= 3){
      res.pb(vi);
      vi.clear();
      vb.pop_back();
      fprintf(stderr,"\n");
    }
  }
  for(int i = 0; i < res.size(); i++){
    res[i] = getpol1(res[i]);
  }
  return res;
}*/


class SmallPolygons{
  
    
  
  
  
  
  
  public:
  vector <string> choosePolygons(vector <int> points, int N){
    np = points.size();
    n = N;
    for(int i = 0; i < np; i+=2){
      ve.pb(PT(points[i],points[i+1],i/2));
    }
    np/=2;
    
    
    
    
    solution res = getsol1();
    
    //solution s(res);
    //cerr << s.getCost() << endl;
    return res.tovs();
  }

};
int main(){
  int Np; 
  vector<int> points;
  int N;
  int v;
  vector<string> ret;
  cin >> Np;
  for (int i=0; i < Np; i++){
    cin >> v;
    points.pb(v);
  }
  cin >> N;
  
  ret = SmallPolygons().choosePolygons(points, N);
  
  
  
  printf("%d\n", ret.size());
  for (int i=0; i < ret.size(); i++){
    printf("%s\n", ret[i].c_str());
    //fprintf(stderr,"VEC: %s\n", ret[i].c_str());
  }
  fflush(stdout);

  return 0;
}
