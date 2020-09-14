   mesh = open('Z:\VGA_Dropbox\Dropbox\VGALabMembers\PI\CVPR2015\Experiments\Guoxian\Data\McGill_noise_new\McGillDBAll_nd002noise_mat\Teddy-bears\b235.mat');
  mesh = mesh.surface;
   coord =[mesh.X mesh.Y mesh.Z];
   triang = mesh.TRIV;
   nopts = length(coord);
   notrg = length(triang);
   
   ofid = fopen('C:\Users\jx7\Desktop\002bearnoise.vtk','w');
fprintf(ofid, '# vtk DataFile Version 3.0\n');
fprintf(ofid,'vtk output\n');
fprintf(ofid,'ASCII\n');
fprintf(ofid,'DATASET POLYDATA\n');
fprintf(ofid,'POINTS %d float\n', nopts);
fprintf(ofid,'%g %g %g\n', coord');
fprintf(ofid,'POLYGONS %d %d\n', notrg, 4*notrg)
fprintf(ofid,'3 %d %d %d\n', triang'-1);
fprintf(ofid,'\n');
fprintf(ofid,'POINT_DATA %d\n', nopts);
fprintf(ofid,'SCALARS distance_from float\n');
fprintf(ofid,'LOOKUP_TABLE default\n');
%fprintf(ofid,'%g\n', distler);
fclose(ofid);