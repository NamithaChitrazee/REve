sap.ui.define([
   'rootui5/eve7/controller/Main.controller',
   'rootui5/eve7/lib/EveManager'
], function(MainController, EveManager) {
   "use strict";

   return MainController.extend("mu2e.controller.Mu2eMain", {

      onInit: function() {
         MainController.prototype.onInit.apply(this, arguments);
         this.mgr.handle.setReceiver(this);
         this.mgr.RegisterController(this);
      },

      onWebsocketMsg: function(handle, msg, offset) {
         this.mgr.onWebsocketMsg(handle, msg, offset);
      },

      // Called once the EveManager has received and built the full element tree.
      onEveManagerInit: function() {
         MainController.prototype.onEveManagerInit.apply(this, arguments);

         var world = this.mgr.childs[0].childs;
         var pthis = this;

         world.forEach(function(item) {
            if (item._typename === "mu2e::GUI") {
               pthis.gui = item;
               // Register the PostStream callback so the server can trigger a UI refresh.
               pthis.mgr.UT_refresh_event_info = function() {
                  pthis.showEventInfo();
               };
               pthis.showEventInfo();
            }
            if (item._typename === "mu2e::EventDisplayManager") {
               pthis.eventMgr = item;
            }
         });
      },

      // Update the browser tab title and the header label with current Run/Subrun/Event.
      showEventInfo: function() {
         if (!this.gui) return;
         var info = "Run: " + this.gui.runid +
                    "  |  Subrun: " + this.gui.subrunid +
                    "  |  Event: " + this.gui.eventid;
         document.title = info;
         var label = this.byId("eventInfoLabel");
         if (label) label.setText(info);
      },

      nextEvent: function() {
         if (!this.eventMgr) return;
         this.mgr.SendMIR("NextEvent()", this.eventMgr.fElementId, "mu2e::EventDisplayManager");
      },

      quitRoot: function() {
         if (!this.eventMgr) return;
         this.mgr.SendMIR("QuitRoot()", this.eventMgr.fElementId, "mu2e::EventDisplayManager");
      }

   });
});
