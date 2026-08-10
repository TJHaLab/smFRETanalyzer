import static org.junit.jupiter.api.Assertions.assertTrue;

import java.awt.GraphicsEnvironment;
import org.junit.jupiter.api.Test;

class HarnessSmokeTest {

    @Test
    void runsHeadless() {
        assertTrue(GraphicsEnvironment.isHeadless());
    }
}
