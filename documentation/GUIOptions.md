# 💡 Mu2e Event Display GUI

The Mu2e Event Display provides an interactive visualization environment based on the ROOT REve framework. The GUI allows users to control what data is displayed, how it is viewed, and how detailed information is retrieved.

## Toggle Menu

The Toggle Menu typically appears as a side panel in the browser, providing quick access to visibility controls and event-specific filtering options. The entire visualization is structurally divided into two primary $\text{REve}$ Scenes to manage the fixed detector elements separately from the changing event data.

The toggle menu includes everything that is displayed in the browser. Our custom products are divided into two primary Scenes:

* $\text{GeomScene}$ (Geometry Scene): This scene contains static, time-independent elements of the detector. This includes the large, unmoving structures of the Mu2e detector, such as the Solenoids, the Target, and the fixed components of the Tracker and Calorimeter. These elements are loaded once and rarely change for a specific data set. Remember to ensure the data-set and the geometry match i.e. MDC2020 data requires MDC2020 geometry.

* $\text{EventScene}$ (Event Data Scene): This scene holds all dynamic, event-dependent data products. This is where the visualization of particle trajectories, reconstructed tracks, hits, and clusters resides. Every time a new event is loaded, the contents of the $\text{EventScene}$ are cleared and redrawn based on the new data products.

* Visibility: You can toggle the overall visibility of major detector components (e.g., turn off the CRV to focus on the Tracker). The same can be done for the event information using the check markers.


## Zooming and Keyboard shortcuts

The visualization scene supports full 3D interaction for detailed inspection of events.

### Mouse Interaction (3D Viewer):

* Left Mouse Drag: Rotates the 3D scene view around the current center point.

* Right Mouse Drag: Pans the scene (moves the camera laterally).

* Scroll Wheel: The mouse can be used to zoom in and out by changing the camera's distance from the scene.

* Keyboard $\text{'Shift'}$ + Left Click: The shift key can be used to recenter the view on the point you click, setting a new pivot point for rotations.

* Keyboard Shortcuts: While many advanced shortcuts exist within the ROOT REve framework, the most commonly used shortcuts are for navigation and resetting the view. Hitting the 'Home' key typically resets the scene back to its predefined initial orientation (often looking down the Z-axis).

## Text Entry

Text entry fields provide critical, powerful filtering capabilities, although the current implementation may be limited in quantity. We are working to advance the options.

Please reach out if you have use-cases.

## MouseOvers

All products have useful "mouseover" labels (tooltips) that link the visual representation directly back to the underlying data.

* Function: Place your mouse on a track, hit, or cluster, and a tooltip will appear, showing information relating to the features of that product.

* Data Provided: The tooltip displays the element's descriptive $\text{Title}$ string, which is custom-built for each data type.

This feature is essential for debugging, event checking, and quantitative analysis within the display.
