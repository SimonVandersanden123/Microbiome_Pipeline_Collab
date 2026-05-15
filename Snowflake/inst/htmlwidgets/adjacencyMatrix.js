HTMLWidgets.widget({
  name: "adjacencyMatrix",
  type: "output",
  factory: function (el, width, height) {

    // Styling
    d3.selectAll('rect.selection').attr('stroke', '#39000E').attr('fill', '#39000E').attr('fill-opacity', .3);

    // --- Dimensions and SVG
    const margin = { top: 0, right: 10, bottom: 10, left: 60 };
    const padding = 30;
    const labelHeight = 100;
    const barcodeDim = {
      height: 50,
      width: width - margin.left - margin.right,
    };
    const adjacencyDim = {
      height: height - margin.top - padding - barcodeDim.height - margin.bottom,
      width: width - margin.left - margin.right,
    };

    const svgAdjacency = d3
      .select(el)
      .append("svg")
      .attr("width", adjacencyDim.width)
      .attr("height", adjacencyDim.height);
    const svgBarcode = d3
      .select(el)
      .append("svg")
      .attr("width", barcodeDim.width)
      .attr("height", barcodeDim.height);

    const gBarcode = svgBarcode
      .append("g")
      .attr("class", "barcode-group-element")
      .attr("transform", `translate(${margin.left},${0})`);
    const gAdjacency = svgAdjacency
      .append("g")
      .attr("class", "adjacency-group-element")
      .attr("transform", `translate(${margin.left},${margin.top + labelHeight})`);

    const gLabels = svgAdjacency
      .append("g")
      .attr("class", "label-group-element")
      .attr("class", "labels")
      .attr("transform", `translate(${margin.left},${margin.top + labelHeight})`)
      .style("font-family", "Verdana, sans-serif");

    // --- Construct scales
    const xScale = d3.scaleBand().range([0, adjacencyDim.width]);
    const xScaleBar = d3.scaleBand().range([0, adjacencyDim.width]);
    // const yScale = d3.scaleBand().range([0, adjacencyDim.height - padding]);
    const yScale = d3.scaleBand().range([0, adjacencyDim.height - labelHeight]);


    // --- Specify brushes
    const brushX = d3.brushX();
    const brush2D = d3.brush();

    // --- Construct groups
    let cells = gAdjacency
      .append("g")
      .attr("class", "matrix")
      .selectAll("rect");

    // Labels

    let rLabels = gLabels
      .append("g")
      .attr("class", "labels-row")
      .selectAll("text");

    let cLabels = gLabels
      .append("g")
      .attr("class", "labels-col")
      .selectAll("text");

    let tLabel = gLabels
      .append("g")
      .append("text")
      .attr("x", -5)
      .attr("y", -5)
      .style("alignment-baseline", "baseline")
      .style("text-anchor", "end")
      .style("fill", "black")
      .text("ASV →");

    // Barcode
    let column_brush = false;
    let barcode = gBarcode.append("g").selectAll("rect");
    // let barcodeLabels = gBarcode.append("g").selectAll("rect");

    // --- Functions
    // scales
    function sourceScale(_key) {
      if (_key === "alphaObserved")
        return d3
          .scaleSequential(d3.interpolateReds)
          .domain([0, d3.max(rows, (d) => d.alphaObserved)]);
      if (_key === "alphaShannon")
        return d3
          .scaleSequential(d3.interpolateReds)
          .domain([0, d3.max(rows, (d) => d.alphaShannon)]);
      if (_key === "alphaChao1")
        return d3
          .scaleSequential(d3.interpolateReds)
          .domain([0, d3.max(rows, (d) => d.alphaChao1)]);
      if (_key === "alphaInvSimpson")
        return d3
          .scaleSequential(d3.interpolateReds)
          .domain([0, d3.max(rows, (d) => d.alphaInvSimpson)]);
      if (_key === "group")
        return d3
          .scaleOrdinal()
          .domain(groups)
          .range([
            "#1b9e77",
            "#d95f02",
            "#7570b3",
            "#e7298a",
            "#66a61e",
            "#e6ab02",
            "#a6761d",
            "#666666",
          ]);
      if (_key === "beta") return cielab;
      if (_key === "abundance")
        return d3.scaleSequential(d3.interpolateReds).domain([0, 10]);
      return (_) => "black";
    }

    function targetScale(_key) {
      if (_key === "observedIn")
        return d3
          .scaleQuantile()
          .domain(d3.range(10))
          .range(["#ffffcc", "#a1dab4", "#41b6c4", "#2c7fb8", "#253494"]);
      if (_key === "abundance")
        return d3.scaleSequentialPow(d3.interpolateTurbo).domain([0, 0.1]);
      if (_key === "kingdom")
        return d3
          .scaleOrdinal()
          .domain(
            [...new Set(columns.map((d) => d[colorTargetNodes]))].filter(
              (k) => k !== undefined
            )
          )
          .range(d3.schemeTableau10)
          .unknown("gray");
      if (_key === "phylum")
        return d3
          .scaleOrdinal()
          .domain(
            [...new Set(columns.map((d) => d[colorTargetNodes]))].filter(
              (k) => k !== undefined
            )
          )
          .range(d3.schemeTableau10)
          .unknown("gray");
      if (_key === "class")
        return d3
          .scaleOrdinal()
          .domain(
            [...new Set(columns.map((d) => d[colorTargetNodes]))].filter(
              (k) => k !== undefined
            )
          )
          .range(d3.schemeTableau10)
          .unknown("gray");
      return (_) => "black";
    }

    function cellScale(_key) {
      if (_key === "relative")
        return d3.scaleSequential(d3.interpolateReds).domain([0, 10]);
      if (_key === "absolute")
        return d3
          .scaleSequential(d3.interpolateReds)
          .domain([0, d3.max(columns, (d) => d.absolute)]);
      return (_) => "black";
    }

    function undef(x) {
      return x === undefined ? 9999 : x;
    }
    // drawing
    function draw() {
      // All width-height related attributes to update on resize
      cells
        .attr("x", (d) => undef(xScale(d.target)))
        .attr("y", (d) => undef(yScale(d.source)))
        .attr("width", xScale.bandwidth())
        .attr("height", yScale.bandwidth());

      rLabels
        .attr("y", (d) => yScale(d.id) + yScale.bandwidth() / 2)
        .style("font-size", yScale.bandwidth() / 3);

      cLabels
        .attr("y", (d) => undef(xScale(d.id)) + (xScale.bandwidth()/2))
        .style("font-size", Math.min(xScale.bandwidth(), 9));

      tLabel.style("font-size", yScale.bandwidth() / 4);
    }

    function drawBarcode() {
      barcode
        .attr("x", (d) => xScaleBar(d.target))
        .attr("width", xScale.bandwidth());
      // barcodeLabels
      //   .attr("x", (d) => xScaleBar(d.id))
      //   .attr("width", xScale.bandwidth());
    }

    function between(x, min, max, range) {
      return x + range > min && x < max;
    }

    return {
      renderValue: function (x) {
        // --- Parse Data
        let {
          graph: { nodes, links },
          settings,
        } = x;

        const rows = nodes.filter((k) => k.type === "source");
        const columns = nodes
          .filter((k) => k.type === "target")
          .filter((d) => d.id);
        if (rows[0]?.seriation) rows.sort((a, b) => a.seriation - b.seriation);
        if (columns[0]?.seriation)
          columns.sort((a, b) => a.seriation - b.seriation);

        xScale.domain(columns.map((d) => d.id));
        xScaleBar.domain(columns.map((d) => d.id));
        yScale.domain(rows.map((d) => d.id));

        let selectedSamples = [];

        // --- Parse Settings
        const {
          source_color,
          target_color,
          cellColor,
          columnBrush,
          crosstalkGroup,
          crosstalkKey,
        } = settings;
        column_brush = columnBrush

        if (!columnBrush) {
          adjacencyDim.height = height - margin.top - margin.bottom - padding;
          svgBarcode.attr("height", 0);
          svgAdjacency.attr("height", adjacencyDim.height);
          yScale.range([0, adjacencyDim.height - padding]);
        }

        // --- Populate Matrix
        const cellColorScale = cellScale(cellColor);
        cells = cells
          .data(links.filter((l) => l.target))
          .join("rect")
          .attr("fill", (d) => cellColorScale(d[cellColor]));
        const rLabelColorScale = sourceScale(source_color);
        rLabels = rLabels
          .data(rows)
          .join("text")
          .attr("x", -5)
          .style("alignment-baseline", "middle")
          .style("text-anchor", "end")
          .style("fill", (d) => rLabelColorScale(d[source_color]))
          .text((d) => d.id)
          .on("mouseover", function () {
            d3.select(this).attr("x", -10);
          })
          .on("mouseout", function () {
            d3.select(this).attr("x", -5);
          })
          .on("click", function (event, datum) {
            const item = d3.select(this);
            if (!selectedSamples.includes(datum.id)) {
              item.style("font-weight", "bold");
              selectedSamples.push(datum.id);
            } else {
              item.style("font-weight", "normal");
              selectedSamples.splice(
                selectedSamples.findIndex((el) => el === datum.id),
                1
              );
            }
            if (selectedSamples.length == 0) {
              cells.attr("opacity", 1);
            } else {
              cells.attr("opacity", 1);
              cells
                .filter((k) => !selectedSamples.includes(k.source))
                .attr("opacity", 0.3);
            }
            update_selection();
          });
        const targetLabelColorScale = targetScale(target_color);
        cLabels = cLabels
          .data(columns)
          .join("text")
          .attr("x", 5)
          .attr("transform", "rotate(-90)")
          .style("alignment-baseline", "middle")
          .style("text-anchor", "start")
          .style("fill", (d) => targetLabelColorScale(d[target_color]))
          .text((d) =>
            d.species
              ? d.species
              : d.genus
              ? d.genus
              : d.families
              ? d.families
              : d.order
              ? d.order
              : d.class
              ? d.class
              : d.phylum
              ? d.phylum
              : d.kingdom
          );
        draw();

        // --- Populate Barcode
        if (columnBrush) {
          barcode = barcode
            .data(links)
            .join("rect")
            // .attr("y", "30%")
            // .attr("height", "70%")
            .attr("y", 0)
            .attr("height", "100%")
            .attr("fill", "black")
            .attr("opacity", 0.2);
          // barcodeLabels = barcodeLabels
          //   .data(columns)
          //   .join("rect")
          //   .attr("y", 0)
          //   .attr("height", "25%")
          //   .attr("fill", "black");
          drawBarcode();
        }

        // --- Configure Horizontal Brush
        brushX.on("brush end", ({ selection }) => {
          if (!selection) {
            xScale.domain(columns.map((d) => d.id));
            draw();
            return;
          }
          const eachBand = xScaleBar.step();
          const [x0, x1] = selection;
          const i0 = Math.ceil(x0 / eachBand);
          const i1 = Math.floor(x1 / eachBand);
          xScale.domain(columns.slice(i0, i1 + 1).map((d) => d.id));
          draw();
        });
        gBarcode.call(brushX);

        // --- Configure crosstalk
        let brush_selection = null;
        const ct_filter = new crosstalk.FilterHandle();

        ct_filter.setGroup(crosstalkGroup);
        brush2D.on("brush end", ({ selection }) => {
          brush_selection = selection;
          update_selection();
        });
        gAdjacency.call(brush2D);

        function update_selection() {
          // Neither selection types active
          if (!brush_selection && !selectedSamples.length) {
            ct_filter.clear();
            return;
          }
          // Only a sample selection active
          keys = null;
          if (!brush_selection) {
            keys = crosstalkKey.filter((k, i) => {
              return selectedSamples.includes(links[i].source);
            });
          }
          // Only a brush selection active
          if (!selectedSamples.length) {
            const [[x0, y0], [x1, y1]] = brush_selection;
            keys = crosstalkKey.filter((k, i) => {
              return (
                between(yScale(links[i].source), y0, y1, yScale.bandwidth()) &&
                between(xScale(links[i].target), x0, x1, xScale.bandwidth())
              );
            });
          }
          if (brush_selection && selectedSamples.length) {
            // Both selections active
            const [[x0, y0], [x1, y1]] = brush_selection;
            keys = crosstalkKey.filter((k, i) => {
              return (
                selectedSamples.includes(links[i].source) &&
                between(yScale(links[i].source), y0, y1, yScale.bandwidth()) &&
                between(xScale(links[i].target), x0, x1, xScale.bandwidth())
              );
            });
          }
          // console.log("Adjacency", keys);
          ct_filter.set(keys);
        }
      },

      resize: function (width, height) {
        const oldInnerWidth = adjacencyDim.width;
        const oldInnerHeight = adjacencyDim.height;
        // --- rescale Barcode
        if (column_brush) {
          barcodeDim.width = width - margin.left - margin.right;
          svgBarcode
            .attr("width", barcodeDim.width)
            .attr("height", barcodeDim.height);
          xScale.range([0, barcodeDim.width]);
          drawBarcode();
        }

        const selectBrushX = d3.brushSelection(gBarcode.node());
        if (selectBrushX) {
          const cw = barcodeDim.width / oldInnerWidth;
          const nx0 = selectBrushX[0] * cw;
          const nx1 = selectBrushX[1] * cw;
          gBarcode.call(brushX.move, [nx0, nx1]);
        }

        // --- rescale Adjacency
        adjacencyDim.width = width - margin.left - margin.right;
        adjacencyDim.height =
          height - margin.top - margin.bottom - barcodeDim.height - padding;
        svgAdjacency
          .attr("width", adjacencyDim.width)
          .attr("height", adjacencyDim.height);
        xScale.range([0, adjacencyDim.width]);
        yScale.range([0, adjacencyDim.height - padding]);
        draw();

        const selectBrush2D = d3.brushSelection(gAdjacency.node());
        if (selectBrush2D) {
          const cw = adjacencyDim.width / oldInnerWidth;
          const ch = adjacencyDim.height / oldInnerHeight;
          const nx = selectBrush2D[0][0] * cw;
          const ny = selectBrush2D[0][1] * ch;
          const nw = selectBrush2D[1][0] * cw;
          const nh = selectBrush2D[1][1] * ch;
          gAdjacency.call(brush2D.move, [
            [nx, ny],
            [nw, nh],
          ]);
        }
      },
    };
  },
});
