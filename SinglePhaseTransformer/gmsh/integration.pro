
Integration {
	{
		Name Integ;
		Case {
			{
				Type Gauss;
   			Case {
    	 				{ GeoElement Point; NumberOfPoints 1; }
    	 				{ GeoElement Line; NumberOfPoints 4; }
    	 				{ GeoElement Triangle; NumberOfPoints 7; }
    	 				{ GeoElement Quadrangle; NumberOfPoints 7; }
    	 				{ GeoElement Tetrahedron; NumberOfPoints 17; }
    	 				{ GeoElement Hexahedron; NumberOfPoints 21; }
    	 				{ GeoElement Pyramid; NumberOfPoints 21; }
    	 				{ GeoElement Prism; NumberOfPoints 21; }
    		}
    	}
		}
	}
}
