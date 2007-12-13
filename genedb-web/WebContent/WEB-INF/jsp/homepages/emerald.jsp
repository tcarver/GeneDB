<%@ include file="/WEB-INF/jsp/topinclude.jspf" %>

<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1">
<link type="text/css" rel="stylesheet" href="<c:url value="/"/>includes/style/emerald.css">
<title>GeneDB Home Page</title>
</head>
<body>
	<div id="wrap">
		<div id="header" align="right">
			<div id="floating-right">
				<img src="<c:url value="/"/>includes/images/fg.jpg"/>
				<img src="<c:url value="/"/>includes/images/he.jpg"/>
				<img src="<c:url value="/"/>includes/images/pl.jpg"/>
				<img src="<c:url value="/"/>includes/images/ve.gif"/>
			</div>
			<h2>Malaria Annotation Workshop - Oct 2007</h2>
			<div id="floating-left">
				<img src="<c:url value="/"/>includes/images/vi.jpg"/>
				<img src="<c:url value="/"/>includes/images/ba.jpg"/>
				<img src="<c:url value="/"/>includes/images/nm.gif"/>
				<img src="<c:url value="/"/>includes/images/bg.jpg"/>
			</div>
		</div>
	    <div id="nav">
	    	<!-- <ul>
	    		<li><a href="#">Home</a></li>
	    		<li><a href="#">Search</a></li>
	    		<li><a href="#">Organisms</a></li>
	    	</ul> -->
	    </div>
	    <div id="sidebar">
	    	<div id="sb">
		    	<%-- <h3>Column 2</h3>
		    	<ul>
		    		<li><a href="#">Link 1</a></li>
		    		<li><a href="#">Link 2</a></li>
		    		<li><a href="#">Link 3</a></li>
		    		<li><a href="#">Link 4</a></li>
		    	</ul> --%>
		    </div>
		    <div id="sb">
		    	<%-- <h3>Column 3</h3>
		    	<ul>
		    		<li><a href="#">Link 1</a></li>
		    		<li><a href="#">Link 2</a></li>
		    		<li><a href="#">Link 3</a></li>
		    		<li><a href="#">Link 4</a></li>
		    	</ul> --%>
		    </div>
	    </div>
	    <div id="main">
	    	<div id="maininner">
		    	<h3>Welcome to the Malaria Re-Annotation (Continuing) Workshop </h3>
		    	<h4>Other websites</h4>
		    	<ul>
		    		<li><a href="http://www.genedb.org/">GeneDB</a></li>
		    		<li><a href="http://www.plasmodb.org/">PlasmoDB</a></li>
		    		<!-- <li><a href="">OPI</a></li>-->
		    	</ul>
		    	
		    	
		    	<h4>Reports</h4>
				<ul>
					<li><a href="./Search/FindCvByName">Browse Cv Terms</a></li>
					<li><a href="./Search/CvTermByCvName?cvName=genedb_products">Product lists</a></li>
					<li><a href="./TopLevelFeatures">Gene lists ordered by location</a></li>
			</div>
	    </div>
	    <div id="footer">
	    	<p></p>
	    </div>
	</div>
</body>
</html>