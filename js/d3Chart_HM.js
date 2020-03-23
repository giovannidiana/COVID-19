// set the dimensions and margins of the graph
var containerWidth = +d3.select('body').style('width').slice(0, -2)

var margin = {top: 10, right: 30, bottom: 60, left: 60},
	width = containerWidth - margin.left - margin.right,
	height = 400 - margin.top - margin.bottom;

// append the svg object to the body of the page

var svg = d3.select("#areaHM")
               .append("svg")
			   .attr("id","svg")
			   .attr("width", width + margin.left + margin.right)
			   .attr("height", height + margin.top + margin.bottom)
			   .append("g")
			   .attr("transform",
					 "translate(" + margin.left + "," + margin.top + ")");

//Read the data
d3.csv("https://raw.githubusercontent.com/giovannidiana/COVID-19/gh-pages/src/AreaListCsv_hm_2203.csv", function(data) {

	var keys0 = d3.keys(data[0]).slice(1,)
	                           .map(function(d) {return d.slice(0,-4);});
	var NC = keys0.length;
	var keys=[];
	for ( var i=0;i<NC/4;i++){
	    keys.push(keys0[i]);
	}

	data.forEach(function(d) {
					d.t = d3.timeParse("%d %b %Y")(d.t); 
					});

	// Deal with select button
		
	currentGroup="Italy";
	var dropDown = d3.select("#dropdown_HM")
	                     .selectAll("option")
	                     .data(keys)
						 .enter()
						 .append('option')
						 .text(function(d) {return d;})
						 .attr("value", function(d) { return d; })
						 .property("selected",function(d) { return d===currentGroup;});

	// get observations
	var observations = data.filter(function(d){
        if(isNaN(+d[currentGroup+"_INF"])){
            return false;
        } else {
		    return +d[currentGroup+"_INF"];
		}
	})

  // Add X axis
		var x = d3.scaleTime()
		          //.domain(d3.extent(data, function(d) { return d3.timeParse("%d %b %Y")(d.t); }))
		          .domain(d3.extent(data, function(d) { return d.t; }))
				  .range([ 0, width ]);
				  
		var xAxis=svg.append("g")
		             .attr("transform", "translate(0," + height + ")")
					 .call(d3.axisBottom(x))
					 //.selectAll("text")
					 //.style("text-anchor", "end")
					 //.attr("transform","rotate(-65)");
  
		
		// Add Y axis
		var y = d3.scaleLinear()
				  .domain([0, d3.max(data, function(d) { return +d[currentGroup+'_qi3'];})])
				  .range([ height, 0]);

		var yAxis=svg.append("g")
                     .call(d3.axisLeft(y));

 // Add a clipPath: everything out of this area won't be drawn.
  
		var clip = svg.append("defs").append("svg:clipPath")
                      .attr("id", "clip")
					  .append("svg:rect")
					  .attr("width", width )
					  .attr("height", height )
					  .attr("x", 0)
					  .attr("y", 0);

	 // Add brushing
		var brush = d3.brush()  // Add the brush feature using the d3.brush function
                      .extent( [ [0,0], [width,height] ] ) // initialise the brush area: start at 0,0 and finishes at width,height: it means I select the whole graph area
                      .on("end", updateChart); // Each

	// Create the scatter variable: where both the circles and the brush take place
		var scatter = svg.append('g')
                         .attr("clip-path", "url(#clip)");


	// Add dots
		scatter.append('g')
			   .selectAll("dot_inf")
			   .data(observations)
			   .enter()
			   .append("circle")
			      .attr("class","inf")
			      .attr("cx", function (d) { return x(d.t); } )
			      .attr("cy", function (d) { return y(d[currentGroup+'_INF']); } )
			      .attr("r", 2)
			      .style("fill","red")


  // add lines
		scatter.append("path")
		       .attr("class","inf")
               .datum(data)
			   .attr("fill", "red")
			   .attr("stroke", "black")
			   .style("opacity", 0.6)
			   .attr("d", d3.area()
							   .x(function(d) { return x(d.t) })
							   .y0(function(d) { return y(d[currentGroup+'_qi1']) })
							   .y1(function(d) { return y(d[currentGroup+'_qi3']) })
							   );

// Add the brushing
		scatter.append("g")
			   .attr("class", "brush")
			   .call(brush);

  // A function that set idleTimeOut to null
		var idleTimeout;
		function idled() { idleTimeout = null; }

  // A function that update the chart for given boundaries
		
		function updateChart() {
				extent = d3.event.selection;

    			// If no selection, back to initial coordinate. Otherwise, update X axis domain
				if(!extent){
						if (!idleTimeout) return idleTimeout = setTimeout(idled, 350); // This allows to wait a little bit
						x.domain(d3.extent(data, function(d) { return d.t; }))
						y.domain([0, d3.max(data, function(d) { return +d[currentGroup+'_qi3'];})])
				} else {
						x.domain([ x.invert(extent[0][0]), x.invert(extent[1][0]) ]);
						y.domain([ y.invert(extent[1][1]), y.invert(extent[0][1]) ]);
						scatter.select(".brush").call(brush.move, null) // This remove the grey brush area as soon as the selection has been done
				}
		  
		  		// Update axis and circle position
				xAxis.transition().duration(500)
					              .call(d3.axisBottom(x))
				yAxis.transition().duration(500).call(d3.axisLeft(y));
				
				scatter.selectAll("circle.inf")
					   .transition().duration(500)
					   .attr("cx", function (d) { return x(d.t); } )
					   .attr("cy", function (d) { return y(d[currentGroup+'_INF']); } );

				scatter.selectAll("path.inf")
						.transition().duration(500)
						.attr("d", d3.area()
										.x(function(d) { return x(d.t) })
										.y0(function(d) { return y(d[currentGroup+'_qi1']) })
										.y1(function(d) { return y(d[currentGroup+'_qi3']) }));
		
		
		};

		function onSelectionChange(currentGroup) {

				x.domain(d3.extent(data, function(d) { return d.t; }))
				y.domain([0, d3.max(data, function(d) { return +d[currentGroup+'_qi3'];})])
				xAxis.transition().duration(500)
					              .call(d3.axisBottom(x))
					              //.selectAll("text")
					              //.style("text-anchor", "end")
					              //.attr("transform","rotate(-65)");
				yAxis.transition().duration(500).call(d3.axisLeft(y));
				
				scatter.selectAll("circle.inf")
					   .transition().duration(500)
					   .attr("cx", function (d) { return x(d.t); } )
					   .attr("cy", function (d) { return y(d[currentGroup+'_INF']); } );
				
				
				scatter.selectAll("path.inf")
						.transition().duration(500)
						.attr("d", d3.area()
										.x(function(d) { return x(d.t) })
										.y0(function(d) { return y(d[currentGroup+'_qi1']) })
										.y1(function(d) { return y(d[currentGroup+'_qi3']) }));
		
		
		};

		

		// Listen to the slider?
		d3.select("#dropdown_HM")
			.on("change", 
				  function(d){
				     currentGroup = this.value;
					 onSelectionChange(currentGroup)
				  });

		
	})
