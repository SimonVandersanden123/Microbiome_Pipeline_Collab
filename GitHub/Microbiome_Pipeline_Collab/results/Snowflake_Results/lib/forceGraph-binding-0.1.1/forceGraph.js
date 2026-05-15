HTMLWidgets.widget({
  name: "forceGraph",
  type: "output",
  factory: function (el, width, height) {
    const legendContainer = d3.select(el).append("div")
      .attr("class", "legend")
      .style("position", "absolute")
      .style("top", 5)
      .style("left", 5);
    const svg = d3.select(el).append("svg")
      .attr("width", width)
      .attr("height", height)
      .attr("viewBox", [-width / 2, -height / 2, width, height])
      .attr("style", "max-width: 100%; height: auto; height: intrinsic;");
    const g = svg.append("g");
    const zoom = d3.zoom().extent([[0, 0],[width, height],]).on("zoom", zoomed);
    svg.call(zoom);

    function zoomed({ transform }) {
      g.attr("transform", transform);
    }

    return {
      renderValue: function (x) {
        // --- Parse and prepare data
        let {
          graph: { nodes, links },
          settings,
        } = x;

        const N = d3.map(nodes, (d) => d.id).map(intern);
        const LS = d3.map(links, ({ source }) => source).map(intern);
        const LT = d3.map(links, ({ target }) => target).map(intern);

        nodes = d3.map(nodes, (d) => ({
          ...d,
          connectedWith: d.connectedWith.split(','),
        }));
        links = d3.map(links, (d, i) => ({
          source: LS[i],
          target: LT[i],
          absolute: d.absolute,
          relative: d.relative,
        }));

        // --- Parse and prepare Settings
        const {
          source_color,
          target_color,
          groupingVariable,
          nodeStroke,
          nodeStrokeWidth,
          nodeOpacity,
          nodeStrokeOpacity,
          highlightAbundance,
          forceOnGroup,
          crosstalkGroup,
          crosstalkKey,
        } = settings;
        let {  
          nodeRadius_source, 
          nodeRadius_target, 
          linkStroke,
          linkStrokeWidth,
          linkStrokeOpacity,
          nodeTitle 
        } = settings;

        let groups;
        if (groupingVariable !== null) groups = [...new Set(d3.map(nodes, (d) => d[groupingVariable]))];
        let tooltipVariables;
        if (nodeTitle == null) {
          nodeTitle = (d) => d.id;
        } else {
          tooltipVariables = nodeTitle;
          nodeTitle = (d) => tooltipVariables.map(v => { if (d[v]) return `${v}: ${d[v]}`}).join('\n');
        }

        nodeRadius_source = Array.isArray(nodeRadius_source) ? constructMapping(nodeRadius_source) : nodeRadius_source;
        nodeRadius_target = Array.isArray(nodeRadius_target) ? constructMapping(nodeRadius_target) : nodeRadius_target;
        linkStroke = Array.isArray(linkStroke) ? constructMapping(linkStroke) : linkStroke;
        linkStrokeWidth = Array.isArray(linkStrokeWidth) ? constructMapping(linkStrokeWidth) : linkStrokeWidth;
        linkStrokeOpacity = Array.isArray(linkStrokeOpacity) ? constructMapping(linkStrokeOpacity) : linkStrokeOpacity;

        let fixedHighlight = false;
        const targetColoringKeys = ["type", "observedIn", "abundance", "kingdom", "phylum", "class", "genus"];

        // --- Construct forces and simulation
        const forceNode = d3.forceManyBody();
        const forceLink = d3.forceLink(links).id(({ index: i }) => N[i]);

        const simulation = d3
          .forceSimulation(nodes)
          .force("link", forceLink)
          .force("charge", forceNode)
          .force("center", d3.forceCenter())
          .force("r", d3.forceRadial((d) => d.type === "source" ? 250 : 200 / d.observedIn).strength((d) => (d.type === "source" ? 0.3 : 0.1)))
          .on("tick", ticked);

        if (forceOnGroup && groups !== undefined) {
          if (groups.length === 2) {
            simulation.force("x", d3.forceX((d) => d.type === "source" ? (d.group === groups[0] ? -250 : 250) : 0).strength((d) => (d.type === "source" ? 0.7 : 0)));
          } else if (groups.length === 3) {
            simulation.force("x", d3.forceX((d) => d.type === "source" ? d.group === groups[0] ? -250 : d.group === groups[2] ? 250 : 0 : 0).strength((d) => (d.type === "source" ? 0.7 : 0)));
            simulation.force("y", d3.forceY((d) => d.type === "source" ? (d.group === groups[1] ? -250 : 250) : 0).strength((d) => (d.type === "source" ? 0.7 : 0)));
          } else if (groups.length === 4) {
            simulation.force('x', d3.forceX((d) => { if (d.type === "source") { if (d.group === groups[0]|| d.group === groups[2]) return -250; else return 250; } }).strength((d) => (d.type === "source" ? 0.7 : 0)));
            simulation.force('y', d3.forceY((d) => { if (d.type === "source") { if (d.group === groups[1]|| d.group === groups[3]) return -250; else return 250; } }).strength((d) => (d.type === "source" ? 0.7 : 0)));
          }
        }

        // --- Drawing

        svg.on('click', resetHighlight);

        // Links
        let link = linkAttrs(g.append("g").selectAll("line").data(links));

        function linkAttrs(link) {
          return link
            .join("line")
            .attr("stroke", (d) => typeof linkStroke === 'object' ? linkStroke.scale(d[linkStroke.varToMap]) : linkStroke)
            .attr("stroke-width", (d) => typeof linkStrokeWidth === 'object' ? linkStrokeWidth.scale(d[linkStrokeWidth.varToMap]) : linkStrokeWidth)
            .attr("stroke-opacity", (d) => typeof linkStrokeOpacity === 'object' ? linkStrokeOpacity.scale(d[linkStrokeOpacity.varToMap]) : linkStrokeOpacity);
        }

        // Nodes
        let node = nodeAttrs(g.append("g").attr("stroke-width", nodeStrokeWidth).selectAll("circle").data(nodes));

        function nodeAttrs(node) {
          node = node
            .join("circle")
            .attr("r", (d) => d.type === "source" ? typeof nodeRadius_source === 'object' ? nodeRadius_source.scale(d[nodeRadius_source.varToMap]) : nodeRadius_source : typeof nodeRadius_target === 'object' ? nodeRadius_target.scale(d[nodeRadius_target.varToMap]) : nodeRadius_target)
            .attr("fill", (d) => setNodeColor(d))
            .attr("stroke", nodeStroke)
            .attr("opacity", nodeOpacity)
            .attr("stroke-opacity", nodeStrokeOpacity)
            .call(drag(simulation))
            .on("mouseenter", d => highlightConnections(event))
            .on("mouseleave", d => highlightAll())
            .on("click", d => fixedHighlightConnections(event));
          node.append("title").text(nodeTitle);
          return node;
        }
        
        // Legend 
        const sourceLegend = sourceScaleLegend(source_color);
        const targetLegend = targetScaleLegend(target_color);
        if (sourceLegend) legendContainer.node().appendChild(sourceLegend);
        if (targetLegend) legendContainer.node().appendChild(targetLegend);


        // --- Configure crosstalk
        const ct_filter = new crosstalk.FilterHandle();
        ct_filter.setGroup(crosstalkGroup);
        const ct_handle = ct_filter.on("change", ({ value, oldValue, sender }) => {
            const filter_links = value ? value.map((id) => links[+id - 1]) : links;
            const filter_nodes = value ? nodes.filter((k) => filter_links.map((d) => d.source.id).concat(filter_links.map((d) => d.target.id)).includes(k.id)) : nodes;

            // Update the simulation
            simulation.stop();
            simulation.nodes(filter_nodes);
            forceLink.links(filter_links);
            simulation.alpha(1);
            simulation.restart();

            // Update drawing
            node = nodeAttrs(node.data(filter_nodes));
            link = linkAttrs(link.data(filter_links));
            ticked();
          }
        );

        // --- Helper functions

        // Node coloring
        function setNodeColor(d) {
          if (d.type === "source") {
            return source_color === "beta" ? sourceScale("beta")(d.beta_mds_1,d.beta_mds_2,d.beta_mds_3) : sourceScale(source_color)(d[source_color]);
          } else return targetColoringKeys.includes(target_color) ? targetScale(target_color)(d[target_color]) : targetScale(target_color)(target_color.includes(d.species) ? d.species : target_color.includes(d.genus) ? d.genus : target_color.includes(d.family) ? d.family : target_color.includes(d.order) ? d.order : target_color.includes(d.class) ? d.class : target_color.includes(d.phylum) ? d.phylum : target_color.includes(d.kingdom) ? d.kingdom : undefined);
        }

        // scales
        function sourceScale(_key) {
          if (_key === "type" || _key === undefined)
            return d3.scaleOrdinal().domain(["source"]).range(["black"]);
          else if (_key === "alphaObserved")
            return d3.scaleSequential(d3.interpolateReds).domain([0, d3.max(nodes, (d) => d.alphaObserved)]);
          else if (_key === "alphaShannon")
            return d3.scaleSequential(d3.interpolateReds).domain([0, d3.max(nodes, (d) => d.alphaShannon)]);
          else if (_key === "alphaChao1")
            return d3.scaleSequential(d3.interpolateReds).domain([0, d3.max(nodes, (d) => d.alphaChao1)]);
          else if (_key === "alphaInvSimpson")
            return d3.scaleSequential(d3.interpolateReds).domain([0, d3.max(nodes, (d) => d.alphaInvSimpson)]);
          else if (_key === "group")
            return d3.scaleOrdinal().domain(groups).range(["#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d","#666666",]);
          else if (_key === "beta") return cielab;
          else if (_key === "abundance")
            return d3.scaleSequentialPow(d3.interpolateReds).domain([0, 10]);
        }

        function targetScale(_key) {
          if (_key === "type" || _key === "undefined")
            return d3.scaleOrdinal().domain(["target"]).range(["steelblue"]);
          else if (_key === "observedIn")
            return d3.scaleQuantile().domain(d3.range(10)).range(["#ffffcc", "#a1dab4", "#41b6c4", "#2c7fb8", "#253494"]);
          else if (_key === "abundance")
            return d3.scaleSequentialPow(d3.interpolateTurbo).domain([0, 0.1]);
          else if (_key === "kingdom")
            return d3.scaleOrdinal().domain([...new Set(nodes.map((d) => d[colorTargetNodes]))].filter((k) => k !== undefined)).range(d3.schemeTableau10).unknown("gray");
          else if (_key === "phylum")
            return d3.scaleOrdinal().domain([...new Set(nodes.map((d) => d[colorTargetNodes]))].filter((k) => k !== undefined)).range(d3.schemeTableau10).unknown("gray");
          else if (_key === "class")
            return d3.scaleOrdinal().domain([...new Set(nodes.map((d) => d[colorTargetNodes]))].filter((k) => k !== undefined)).range(d3.schemeTableau10).unknown("gray");
          else if (_key === "genus") {
  // Create the scale
          const genusDomain = [...new Set(nodes.map((d) => d[colorTargetNodes]))].filter((k) => k !== undefined);

  // Check if the 'Target_Genus' exists in the domain
           const Target_Genus = "Target_Genus"; // Replace this with the actual genus name or dynamic input as needed

  // Create a color scale
           const colorScale = d3.scaleOrdinal()
            .domain(genusDomain)
            .range(d3.schemeTableau10)
            .unknown("gray");

  // Apply red color for 'Target_Genus'
            return (d) => {
            if (d[colorTargetNodes] === Target_Genus) {
             return "red"; // Color Target_Genus in red
            } else {
             return colorScale(d); // Use the default color scale for other genera
            }
           };
}

            return d3.scaleOrdinal().domain(_key).range(d3.schemePaired).unknown("gray");
        }

        function cielab(x, y, z) {
          const domain = [-1, 1];
          const c_l = d3.scaleLinear().domain(domain).range([0, 100]);
          const c_ab = d3.scaleLinear().domain(domain).range([-160, 160]);
          const [l, a, b] = [c_l(z), c_ab(y), c_ab(x)];
          return d3.lab(l, a, b);
        }

        function constructMapping(_arr) {

          let [variable, scaleType, domain, range, unknown] = _arr;

          let mapping;
          switch (scaleType) {
            case "linear":
              if (typeof range[0] === 'number') {
                mapping = d3.scaleLinear().domain(domain).range(range).clamp(true);
              } else {
                mapping = d3.scaleLinear().domain(domain).range(range);
              }
              break;
            case "ordinal":
              mapping = d3.scaleOrdinal().domain(domain).range(range);
              break;
            case "sequential":
              mapping = d3.scaleSequential().domain(domain).range(range);
              break;
            case "quantile":
              mapping = d3.scaleQuantile().domain(domain).range(range);
              break;
          }
          
          return {
            varToMap: variable,
            scale: mapping
          }

        }

        // preparation
        function intern(value) {
          return value !== null && typeof value === "object" ? value.valueOf() : value;
        }

        // Dragging
        function drag(simulation) {
          function dragstarted(event) {
            if (!event.active) simulation.alphaTarget(0.3).restart();
            event.subject.fx = event.subject.x;
            event.subject.fy = event.subject.y;
          }

          function dragged(event) {
            event.subject.fx = event.x;
            event.subject.fy = event.y;
          }

          function dragended(event) {
            if (!event.active) simulation.alphaTarget(0);
            event.subject.fx = null;
            event.subject.fy = null;
          }

          return d3.drag()
            .on("start", dragstarted)
            .on("drag", dragged)
            .on("end", dragended);
        }

        // Drawing
        function ticked() {
          link
            .attr("x1", (d) => d.source.x)
            .attr("y1", (d) => d.source.y)
            .attr("x2", (d) => d.target.x)
            .attr("y2", (d) => d.target.y);

          node.attr("cx", (d) => d.x).attr("cy", (d) => d.y);
        }

        // Node highlighting
        function highlightConnections(event) {
          const selected = event.target.__data__;
          const nodesToHighlight = selected.connectedWith.concat(selected.id);
          const linksToHighlight = links.filter(k => k.target.id === selected.id);
          node.filter(k => !nodesToHighlight.includes(k.id)).attr('opacity', .2);
          if (selected.type === "target" && highlightAbundance) node.filter(k => nodesToHighlight.includes(k.id) && k.type === "source").attr("fill", d => sourceScale("abundance")(linksToHighlight.filter((l) => l.source === d)[0].relative));
          link.filter(k => !(nodesToHighlight.includes(k.source.id) && nodesToHighlight.includes(k.target.id))).attr('opacity', .2);
        }
        function highlightAll() {
          node.attr('opacity', d => d.type === 'source' ? nodeOpacity : nodeOpacity - .1)
            .attr("stroke-opacity", nodeStrokeOpacity)
            .attr('fill', (d) => setNodeColor(d))
            .attr("stroke", nodeStroke);
          link
            .attr('opacity', (d) => typeof linkStrokeOpacity === 'object' ? linkStrokeOpacity.scale(d[linkStrokeOpacity.varToMap]) : linkStrokeOpacity)
            .attr('stroke', (d) => typeof linkStroke === 'object' ? linkStroke.scale(d[linkStroke.varToMap]) : linkStroke);
        }
        function fixedHighlightConnections(event) {
          fixedHighlight = true;
          highlightAll();
          highlightConnections(event);
        }
        function resetHighlight(event) {
          if (event.target === svg.node()) {
            fixedHighlight = false;
            highlightAll();
          }
        }

        function sourceScaleLegend(_key) {
          if (_key === "type" || _key === undefined)
            return undefined;
          else if (_key === "alphaObserved")
            return Legend(d3.scaleSequential(d3.interpolateReds).domain([0, d3.max(nodes, (d) => d.alphaObserved)]), { title: `Source Nodes: ${source_color}` });
          else if (_key === "alphaShannon")
            return Legend(d3.scaleSequential(d3.interpolateReds).domain([0, d3.max(nodes, (d) => d.alphaShannon)]), { title: `Source Nodes: ${source_color}` });
          else if (_key === "alphaChao1")
            return Legend(d3.scaleSequential(d3.interpolateReds).domain([0, d3.max(nodes, (d) => d.alphaChao1)]), { title: `Source Nodes: ${source_color}` });
          else if (_key === "alphaInvSimpson")
            return Legend(d3.scaleSequential(d3.interpolateReds).domain([0, d3.max(nodes, (d) => d.alphaInvSimpson)]), { title: `Source Nodes: ${source_color}` });
          else if (_key === "group")
            return Legend(d3.scaleOrdinal().domain(groups).range(["#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d","#666666",]), { title: `Source Nodes: ${source_color}` });
          // else if (_key === "beta") return cielab;
          else if (_key === "abundance")
            return Legend(d3.scaleSequentialPow(d3.interpolateReds).domain([0, 10]), { title: `Source Nodes: ${source_color}` });
        }

        function targetScaleLegend(_key) {
          if (_key === "type" || _key === "undefined")
            return undefined;
          else if (_key === "observedIn")
            return Legend(d3.scaleQuantile().domain(d3.range(10)).range(["#ffffcc", "#a1dab4", "#41b6c4", "#2c7fb8", "#253494"]), { title: `Target Nodes: ${target_color}` });
          else if (_key === "abundance")
            return Legend(d3.scaleSequentialPow(d3.interpolateTurbo).domain([0, 0.1]), { title: `Target Nodes: ${target_color}` });
          else if (_key === "kingdom")
            return Legend(d3.scaleOrdinal().domain([...new Set(nodes.map((d) => d[colorTargetNodes]))].filter((k) => k !== undefined)).range(d3.schemeTableau10).unknown("gray"), { title: `Target Nodes: ${target_color}` });
          else if (_key === "phylum")
            return Legend(d3.scaleOrdinal().domain([...new Set(nodes.map((d) => d[colorTargetNodes]))].filter((k) => k !== undefined)).range(d3.schemeTableau10).unknown("gray"), { title: `Target Nodes: ${target_color}` });
          else if (_key === "class")
            return Legend(d3.scaleOrdinal().domain([...new Set(nodes.map((d) => d[colorTargetNodes]))].filter((k) => k !== undefined)).range(d3.schemeTableau10).unknown("gray"), { title: `Target Nodes: ${target_color}` });
          else if (_key === "genus")
            return Legend(d3.scaleOrdinal().domain([...new Set(nodes.map((d) => d[colorTargetNodes]))].filter((k) => k !== undefined)).range(d3.schemeTableau10).unknown("gray"), { title: `Target Nodes: ${target_color}` });
          else 
            return Legend(d3.scaleOrdinal().domain(_key).range(d3.schemePaired).unknown("gray"), { title: `Target Nodes: ${target_color}` });
        }

        function Swatches(color, {
          title
        } = {}) {
          
          const svg = d3.create("svg")
              .attr("width", width)
              .attr("height", height)
              .attr("viewBox", [0, 0, width, height])
              .style("background-color", "white")
              .style("overflow", "visible")
              .style("display", "block");
        }

        function Legend(color, {
          title,
          tickSize = 6,
          width = 320, 
          height = 44 + tickSize,
          marginTop = 18,
          marginRight = 0,
          marginBottom = 16 + tickSize,
          marginLeft = 0,
          ticks = width / 64,
          tickFormat,
          tickValues
        } = {}) {
        
          function ramp(color, n = 256) {
            const canvas = document.createElement("canvas");
            canvas.width = n;
            canvas.height = 1;
            const context = canvas.getContext("2d");
            for (let i = 0; i < n; ++i) {
              context.fillStyle = color(i / (n - 1));
              context.fillRect(i, 0, 1, 1);
            }
            return canvas;
          }
        
          const svg = d3.create("svg")
              .attr("width", width)
              .attr("height", height)
              .attr("viewBox", [0, 0, width, height])
              .style("background-color", "white")
              .style("overflow", "visible")
              .style("display", "block");
        
          let tickAdjust = g => g.selectAll(".tick line").attr("y1", marginTop + marginBottom - height);
          let x;
        
          // Continuous
          if (color.interpolate) {
            const n = Math.min(color.domain().length, color.range().length);
        
            x = color.copy().rangeRound(d3.quantize(d3.interpolate(marginLeft, width - marginRight), n));
        
            svg.append("image")
                .attr("x", marginLeft)
                .attr("y", marginTop)
                .attr("width", width - marginLeft - marginRight)
                .attr("height", height - marginTop - marginBottom)
                .attr("preserveAspectRatio", "none")
                .attr("xlink:href", ramp(color.copy().domain(d3.quantize(d3.interpolate(0, 1), n))).toDataURL());
          }
        
          // Sequential
          else if (color.interpolator) {
            x = Object.assign(color.copy()
                .interpolator(d3.interpolateRound(marginLeft, width - marginRight)),
                {range() { return [marginLeft, width - marginRight]; }});
        
            svg.append("image")
                .attr("x", marginLeft)
                .attr("y", marginTop)
                .attr("width", width - marginLeft - marginRight)
                .attr("height", height - marginTop - marginBottom)
                .attr("preserveAspectRatio", "none")
                .attr("xlink:href", ramp(color.interpolator()).toDataURL());
        
            // scaleSequentialQuantile doesn’t implement ticks or tickFormat.
            if (!x.ticks) {
              if (tickValues === undefined) {
                const n = Math.round(ticks + 1);
                tickValues = d3.range(n).map(i => d3.quantile(color.domain(), i / (n - 1)));
              }
              if (typeof tickFormat !== "function") {
                tickFormat = d3.format(tickFormat === undefined ? ",f" : tickFormat);
              }
            }
          }
        
          // Threshold
          else if (color.invertExtent) {
            const thresholds
                = color.thresholds ? color.thresholds() // scaleQuantize
                : color.quantiles ? color.quantiles() // scaleQuantile
                : color.domain(); // scaleThreshold
        
            const thresholdFormat
                = tickFormat === undefined ? d => d
                : typeof tickFormat === "string" ? d3.format(tickFormat)
                : tickFormat;
        
            x = d3.scaleLinear()
                .domain([-1, color.range().length - 1])
                .rangeRound([marginLeft, width - marginRight]);
        
            svg.append("g")
              .selectAll("rect")
              .data(color.range())
              .join("rect")
                .attr("x", (d, i) => x(i - 1))
                .attr("y", marginTop)
                .attr("width", (d, i) => x(i) - x(i - 1))
                .attr("height", height - marginTop - marginBottom)
                .attr("fill", d => d);
        
            tickValues = d3.range(thresholds.length);
            tickFormat = i => thresholdFormat(thresholds[i], i);
          }
        
          // Ordinal
          else {
            x = d3.scaleBand()
                .domain(color.domain())
                .rangeRound([marginLeft, width - marginRight]);
        
            svg.append("g")
              .selectAll("rect")
              .data(color.domain())
              .join("rect")
                .attr("x", x)
                .attr("y", marginTop)
                .attr("width", Math.max(0, x.bandwidth() - 1))
                .attr("height", height - marginTop - marginBottom)
                .attr("fill", color);
        
            tickAdjust = () => {};
          }
        
          svg.append("g")
              .attr("transform", `translate(0,${height - marginBottom})`)
              .call(d3.axisBottom(x)
                .ticks(ticks, typeof tickFormat === "string" ? tickFormat : undefined)
                .tickFormat(typeof tickFormat === "function" ? tickFormat : undefined)
                .tickSize(tickSize)
                .tickValues(tickValues))
              .call(tickAdjust)
              .call(g => g.select(".domain").remove())
              .call(g => g.append("text")
                .attr("x", marginLeft)
                .attr("y", marginTop + marginBottom - height - 6)
                .attr("fill", "currentColor")
                .attr("text-anchor", "start")
                .attr("font-weight", "bold")
                .attr("class", "title")
                .text(title));
        
          return svg.node();
        }
      },

      resize: function (width, height) {
        svg.attr("width", width).attr("height", height)
          .attr("viewBox", [-width / 2, -height / 2, width, height]);

        zoom.extent([[0, 0],[width, height],]);
      },
    };
  },
});
