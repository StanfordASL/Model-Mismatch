function teb_box = define_teb_obj(teb)

teb_box_V = [-teb;
           -teb+2*[teb(1),0,0];
           -teb+2*[teb(1),teb(2),0];
           -teb+2*[0,teb(2),0];
           -teb+2*[0,0,teb(3)];
           -teb+2*[teb(1),0,teb(3)];
           -teb+2*[teb(1),teb(2),teb(3)];
           -teb+2*[0,teb(2),teb(3)]];           

box = convhull(teb_box_V(:,1),teb_box_V(:,2),teb_box_V(:,3));

teb_box = struct('hull',box,'V',teb_box_V,'nV',size(teb_box_V,1));

end