HTMLWidgets.widget({
    name: "barcode",
    type: "output",
    factory: function (el, width, height) {
      // --- Dimensions and SVG
      const margin = { top: 10, right: 10, bottom: 10, left: 10 };
      let innerHeight = height - margin.top - margin.bottom;
      let innerWidth = width - margin.left - margin.right;
  
      const svg = d3
        .select(el)
        .append("svg")
        .attr("width", innerWidth)
        .attr("height", innerHeight);
      
      const g = svg
        .append("g")
        .attr("class", "barcode-group-element")
        .attr("transform", `translate(${margin.left},${margin.top})`);
  
  
      // --- Construct scales
      const xScale = d3.scaleBand().range([0, innerWidth]);
    //   const colorScale = 
  
      // --- Specify brush
      const brushX = d3.brushX();
  
      // Barcode
      let barcode = g.append("g")
        .selectAll("rect");
      let barcodeLabels = g.append("g")
        .selectAll("rect");
  
      function drawBarcode() {
        barcode
          .attr("x", d => xScale(d.target))
          .attr("width", xScale.bandwidth());
        barcodeLabels
          .attr("x", d => xScale(d.id))
          .attr("width", xScale.bandwidth());
      }
  
      return {
        renderValue: function (x) {
          // --- Parse Data
          let {
            graph: { nodes, links },
            settings,
          } = x;
  
          const columns = nodes.filter((k) => k.type === "target");
          if (columns[1].seriation) columns.sort((a,b) => a.seriation - b.seriation)
  
          xScale.domain(columns.map((d) => d.id));
  
          // --- Parse Settings
          const {
            color,
            crosstalkGroup,
            crosstalkKey,
          } = settings;
  
          // --- Populate Barcode
          barcode = barcode
            .data(links)
            .join("rect")
            .attr("y", "30%")
            .attr("height", "75%")
            .attr("fill", color ? color : "black")
            .attr("opacity", 0.2);
          barcodeLabels = barcodeLabels
            .data(columns)
            .join("rect")
            .attr("y", 0)
            .attr("height", "25%")
            .attr("fill", "black");
          drawBarcode();
  
          // --- Configure Horizontal Brush
          const ct_filter = new crosstalk.FilterHandle();
          ct_filter.setGroup(crosstalkGroup);
          brushX.on("brush end", ({selection}) => {
            if(selection) {
              const [x0, x1] = selection;
              let selectedColumns = columns.filter(k => x0 <= xScale(k.id) && xScale(k.id) < x1 || x0 <= (xScale(k.id) + xScale.bandwidth()) && xScale(k.id) < x1).map(d => d.id);
              const keys = crosstalkKey.filter(
                (k,i) => {
                    selectedColumns.includes(links[i].target)
                }
              );
              ct_filter.set(keys);
            } else {
                ct_filter.clear();
            }
          });
          g.call(brushX)
  
        },
  
        resize: function (width, height) {
  
          const oldInnerWidth = innerWidth;
          const oldInnerHeight = innerHeight;
          // --- rescale Barcode
          svg.attr("width", width).attr("height", height);
          innerWidth = width - margin.left - margin.right;
          xScale.range([0, innerWidth]);
          drawBarcode();
        
          const select = d3.brushSelection(g.node());
          if (select) {
            const cw = innerWidth / oldInnerWidth;
            const nx0 = select[0] * cw;
            const nx1 = select[1] * cw;
            g.call(brushX.move, [nx0,nx1]);
          }
        },
      };
    },
  });
  