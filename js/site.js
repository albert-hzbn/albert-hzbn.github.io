// Shared page behaviour: scroll reveals, header shadow, figure lightbox and the
// screenshot showcase. Everything degrades to static content without JS and
// honours prefers-reduced-motion.
(function () {
    const reduceMotion = window.matchMedia('(prefers-reduced-motion: reduce)').matches;

    // ---------- Header shadow once the page scrolls ----------
    const header = document.querySelector('.site-header');
    if (header) {
        const onScroll = () => header.classList.toggle('scrolled', window.scrollY > 8);
        window.addEventListener('scroll', onScroll, { passive: true });
        onScroll();
    }

    // ---------- Scroll reveals ----------
    const revealables = document.querySelectorAll('[data-reveal]');
    if (reduceMotion || !('IntersectionObserver' in window)) {
        revealables.forEach((el) => el.classList.add('in'));
    } else {
        // Children of a [data-reveal-stagger] parent get increasing delays.
        document.querySelectorAll('[data-reveal-stagger]').forEach((parent) => {
            Array.from(parent.children).forEach((child, i) => {
                child.setAttribute('data-reveal', '');
                child.style.setProperty('--reveal-delay', (i * 70) + 'ms');
            });
        });
        const io = new IntersectionObserver((entries) => {
            entries.forEach((entry) => {
                if (entry.isIntersecting) {
                    entry.target.classList.add('in');
                    io.unobserve(entry.target);
                }
            });
        }, { rootMargin: '0px 0px -8% 0px', threshold: 0.08 });
        document.querySelectorAll('[data-reveal]').forEach((el) => io.observe(el));
    }

    // ---------- Figure lightbox ----------
    const zoomables = document.querySelectorAll('.figure img[data-zoom]');
    if (zoomables.length) {
        const box = document.createElement('div');
        box.className = 'lightbox';
        box.setAttribute('role', 'dialog');
        box.setAttribute('aria-modal', 'true');
        box.setAttribute('aria-label', 'Enlarged figure');
        box.innerHTML = '<figure><img alt=""><figcaption></figcaption></figure>' +
            '<button class="lightbox-close" aria-label="Close">×</button>';
        document.body.appendChild(box);
        const bigImg = box.querySelector('img');
        const bigCap = box.querySelector('figcaption');
        let lastFocus = null;

        const open = (img) => {
            lastFocus = document.activeElement;
            bigImg.src = img.currentSrc || img.src;
            bigImg.alt = img.alt;
            const cap = img.closest('figure').querySelector('figcaption');
            bigCap.textContent = cap ? cap.textContent : '';
            box.classList.add('open');
            box.querySelector('.lightbox-close').focus();
        };
        const close = () => {
            box.classList.remove('open');
            if (lastFocus) lastFocus.focus();
        };

        zoomables.forEach((img) => {
            img.tabIndex = 0;
            img.setAttribute('role', 'button');
            img.setAttribute('aria-label', 'Enlarge figure: ' + img.alt);
            img.addEventListener('click', () => open(img));
            img.addEventListener('keydown', (e) => {
                if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); open(img); }
            });
        });
        box.addEventListener('click', (e) => { if (e.target !== bigImg) close(); });
        document.addEventListener('keydown', (e) => {
            if (e.key === 'Escape' && box.classList.contains('open')) close();
        });
    }

    // ---------- Screenshot showcase (cross-fade) ----------
    document.querySelectorAll('.showcase').forEach((show) => {
        const slides = Array.from(show.querySelectorAll('.slide'));
        const caption = show.querySelector('.showcase-caption');
        const dots = show.querySelector('.showcase-dots');
        const progress = show.querySelector('.showcase-progress');
        if (slides.length < 2) return;
        let current = 0, timer = null;

        const buttons = slides.map((slide, i) => {
            const b = document.createElement('button');
            b.type = 'button';
            b.setAttribute('aria-label', 'Show screenshot ' + (i + 1) + ': ' + slide.dataset.caption);
            b.addEventListener('click', () => { go(i); restart(); });
            dots.appendChild(b);
            return b;
        });

        function go(i) {
            slides[current].classList.remove('active');
            buttons[current].removeAttribute('aria-current');
            current = (i + slides.length) % slides.length;
            slides[current].classList.add('active');
            buttons[current].setAttribute('aria-current', 'true');
            caption.textContent = slides[current].dataset.caption;
            if (progress) {
                progress.classList.remove('run');
                void progress.offsetWidth; // restart the CSS animation
                if (!reduceMotion) progress.classList.add('run');
            }
        }
        function restart() {
            clearInterval(timer);
            if (!reduceMotion) timer = setInterval(() => go(current + 1), 5000);
        }

        go(0);
        restart();
        show.addEventListener('mouseenter', () => { clearInterval(timer); show.classList.add('paused'); });
        show.addEventListener('mouseleave', () => { show.classList.remove('paused'); go(current); restart(); });
    });
})();
