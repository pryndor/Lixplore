
document.addEventListener('DOMContentLoaded', () => {
    document.querySelectorAll('.toggle-btn').forEach(btn => {
        btn.addEventListener('click', function() {
            const card = this.closest('.result-card'); // target the correct card
            card.classList.toggle('expanded');

            // Optionally, change button text
            if (card.classList.contains('expanded')) {
                this.textContent = 'Show less';
            } else {
                this.textContent = 'Show more';
            }
        });
    });
});




