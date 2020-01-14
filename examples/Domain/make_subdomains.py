import SVMTK as svm


if __name__ == "__main__":
    
 
   sf= svm.SubdomainMap() 
   s = []
   for i in range(4):
      print(i)
      s.append(svm.Surface())
      s[-1].make_sphere(0,0,0,25.-3.*i)
      print(s[-1])



   sf.add("1000",1)
   sf.add("1100",2)  
   sf.add("1110",3) 
   sf.add("1111",4) 
   sf.print()
   maker = svm.Domain(s,sf)

   
   maker.create_mesh(32)

   maker.save("subdomains.mesh")
   boundary = maker.get_boundary(0)
   boundary.save("subdomain_boundary.off")
