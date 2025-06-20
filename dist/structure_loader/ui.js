// ui.js
document.addEventListener("DOMContentLoaded", () => {
    const collapsibles = document.querySelectorAll(".collapsible");
  
    collapsibles.forEach(button => {
      button.addEventListener("click", function () {
        this.classList.toggle("active");
        const content = this.nextElementSibling;
  
        if (content.style.maxHeight) {
          content.style.maxHeight = null;
        } else {
          content.style.maxHeight = content.scrollHeight + "px";
        }
      });
    });
  });
  