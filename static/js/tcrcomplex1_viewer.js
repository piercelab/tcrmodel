// Create NGL Stage object
var stage = new NGL.Stage( "nglviewport" );
stage.setParameters({ backgroundColor: "white" });

// Handle window resizing
window.addEventListener( "resize", function( event ){
    stage.handleResize();
}, false );


// Code for example: interactive/simple-viewer
function addElement (el) {
  Object.assign(el.style, {
    position: "absolute",
    zIndex: 10
  })
  stage.viewer.container.appendChild(el)
}

function createElement (name, properties, style) {
  var el = document.createElement(name)
  Object.assign(el, properties)
  Object.assign(el.style, style)
  return el
}

function createSelect (options, properties, style) {
  var select = createElement("select", properties, style)
  options.forEach(function (d) {
    select.add(createElement("option", {
      value: d[ 0 ], text: d[ 1 ]
    }))
  })
  return select
}

function loadStructure (input) {
  stage.removeAllComponents()
  return stage.loadFile(input).then(function (o) {
    o.autoView()
    o.addRepresentation(polymerSelect.value, {
      sele: "polymer",
      name: "polymer"
    })
    o.addRepresentation("ball+stick", {
      name: "ligand",
      visible: ligandCheckbox.checked,
      sele: "not ( polymer or water or ion )"
    })
    o.addRepresentation("spacefill", {
      name: "waterIon",
      visible: waterIonCheckbox.checked,
      sele: "water or ion",
      scale: 0.25
    })
  })
}

var polymerSelect = createSelect([
  [ "cartoon", "cartoon" ],
  [ "spacefill", "spacefill" ],
  [ "licorice", "licorice" ],
  [ "surface", "surface" ]
], {
  onchange: function (e) {
    stage.getRepresentationsByName("polymer").dispose()
    stage.eachComponent(function (o) {
      o.addRepresentation(e.target.value, {
        sele: "polymer",
        name: "polymer"
      })
    })
  }
}, { top: "36px", left: "12px" })
//addElement(polymerSelect)

var centerButton = createElement("input", {
  type: "button",
  value: "center",
  onclick: function () {
    stage.autoView(1000)
  }
}, { top: "108px", left: "12px" })
//addElement(centerButton)

function nglCenter() {
    stage.autoView(1000);
}

function polymerRender (e) {
    stage.getRepresentationsByName("polymer").dispose()
    stage.eachComponent(function (o) {
      o.addRepresentation(e, {
        sele: "polymer",
        name: "polymer"
      })
    })
  }

$("#ddsec a").on('click', function(e) {
  e.preventDefault(); // cancel the link behaviour
  var selText = $(this).text();
    polymerRender(selText);
});

loadStructure(modelfname)

